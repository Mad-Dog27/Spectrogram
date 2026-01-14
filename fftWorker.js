/* fftWorker.js
   Purpose:
   - Receive resampled audio stream (Float32 chunks) from micWorker
   - Maintain FIFO
   - Produce STFT frames with overlap via hop size:
       hop = frameSize - round(frameSize * overlapPercent)
   - Apply chosen window
   - Zero-pad to NFFT
   - FFT
   - Output magnitude or dB (EPS-safe, normalized, clamped)

   Output:
   - Float32Array length NFFT/2 (positive-frequency bins)
*/

"use strict";

let SAMPLE_RATE = 16000;
let FRAME_SIZE = 256;
let OVERLAP_PERCENT = 0.75;
let NFFT = 2048;

let CHOSEN_WINDOW = "hamming"; // "rectangular" | "hamming" | "blackmanHarris"
let MAG_MODE = "magnitude";    // "magnitude" | "deciBels"

let PAUSED = false;

// dB display settings
let DB_REF = 1;       // reference AFTER normalization (keep 1 if using normalized magnitude)
let DB_MIN = -40;     // clamp range for display
let DB_MAX = 0;

const EPS = 1e-12;

// stream FIFO at target rate
let fifo = new Float32Array(0);

function appendFloat32(a, b) {
  if (a.length === 0) return b;
  if (b.length === 0) return a;
  const out = new Float32Array(a.length + b.length);
  out.set(a, 0);
  out.set(b, a.length);
  return out;
}

function getHop() {
  const overlap = Math.round(FRAME_SIZE * OVERLAP_PERCENT);
  return Math.max(1, FRAME_SIZE - overlap);
}

function applyWindow(frame) {
  const N = frame.length;
  const out = new Float32Array(N);

  if (CHOSEN_WINDOW === "rectangular") {
    out.set(frame);
    return out;
  }

  if (CHOSEN_WINDOW === "hamming") {
    for (let n = 0; n < N; n++) {
      const w = 0.54 - 0.46 * Math.cos((2 * Math.PI * n) / (N - 1));
      out[n] = frame[n] * w;
    }
    return out;
  }

  // blackmanHarris
  for (let n = 0; n < N; n++) {
    const a0 = 0.35875;
    const a1 = 0.48829;
    const a2 = 0.14128;
    const a3 = 0.01168;
    const w =
      a0
      - a1 * Math.cos((2 * Math.PI * n) / (N - 1))
      + a2 * Math.cos((4 * Math.PI * n) / (N - 1))
      - a3 * Math.cos((6 * Math.PI * n) / (N - 1));
    out[n] = frame[n] * w;
  }
  return out;
}

function zeroPad(frame) {
  const out = new Float32Array(NFFT);
  // place at start (simple). Centering not needed for magnitude spectra.
  out.set(frame.slice(0, Math.min(frame.length, NFFT)), 0);
  return out;
}

function bitReverseIndex(index, bits) {
  let reversed = 0;
  for (let i = 0; i < bits; i++) {
    reversed = (reversed << 1) | (index & 1);
    index >>= 1;
  }
  return reversed;
}

function fftReal(input) {
  const N = input.length;
  const levels = Math.log2(N);
  if ((N & (N - 1)) !== 0) throw new Error("FFT length must be power of 2");

  const out = new Array(N);
  for (let i = 0; i < N; i++) {
    const j = bitReverseIndex(i, levels);
    out[i] = { real: input[j], imag: 0 };
  }

  for (let size = 2; size <= N; size <<= 1) {
    const half = size >> 1;
    const step = (-2 * Math.PI) / size;

    for (let i = 0; i < N; i += size) {
      for (let j = 0; j < half; j++) {
        const even = out[i + j];
        const odd = out[i + j + half];

        const ang = step * j;
        const wr = Math.cos(ang);
        const wi = Math.sin(ang);

        const tr = wr * odd.real - wi * odd.imag;
        const ti = wr * odd.imag + wi * odd.real;

        out[i + j] = { real: even.real + tr, imag: even.imag + ti };
        out[i + j + half] = { real: even.real - tr, imag: even.imag - ti };
      }
    }
  }

  return out;
}

// Normalized single-sided magnitude (amplitude-ish)
function magSingleSided(fftBins) {
  const N = fftBins.length;
  const half = N >> 1;
  const mags = new Float32Array(half);

  for (let i = 0; i < half; i++) {
    let mag = Math.hypot(fftBins[i].real, fftBins[i].imag);

    // Normalize by NFFT (keeps magnitude stable across NFFT)
    mag /= N;

    // Single-sided correction (except DC and Nyquist)
    if (i !== 0 && i !== half) mag *= 2;

    mags[i] = mag;
  }
  return mags;
}

let refEMA = 1e-4; // start at MIN_REF

function toDB(mags) {
  const out = new Float32Array(mags.length);
  const MIN_REF = 1e-4;
  const DB_MIN = -80;

  // 95th percentile reference
  const sorted = Array.from(mags).sort((a,b)=>a-b);
  let frameRef = sorted[Math.floor(0.95 * sorted.length)];
  frameRef = Math.max(frameRef, MIN_REF);

  // smooth reference (EMA)
  refEMA = Math.max(
    MIN_REF,
    0.15 * frameRef + 0.85 * refEMA
  );

  for (let i = 0; i < mags.length; i++) {
    const db = 20 * Math.log10((mags[i] + 1e-12) / refEMA);
    out[i] = Math.max(DB_MIN, Math.min(0, db));
  }
  return out;
}



function processFrames() {
  const hop = getHop();


    
  while (fifo.length >= FRAME_SIZE) {
    const frame = fifo.slice(0, FRAME_SIZE);
    fifo = fifo.slice(hop);

    const windowed = applyWindow(frame);
    const padded = zeroPad(windowed);
    const bins = fftReal(padded);
    const mags = magSingleSided(bins);

    let out;
    if (MAG_MODE === "deciBels") out = toDB(mags);
    else out = mags;
    if (MAG_MODE === "deciBels") {
  // after mags computed
  let mn = 1e9, mx = -1e9;
  for (let i=0;i<mags.length;i++){ const v=mags[i]; if(v<mn) mn=v; if(v>mx) mx=v; }
  postMessage({ type:"debug", magsMin: mn, magsMax: mx, dbRef: DB_REF });
}
    postMessage(out, [out.buffer]);
  }
}

onmessage = (e) => {
  const msg = e.data;

  if (msg && msg.type === "config") {
    if (typeof msg.sampleRate === "number") SAMPLE_RATE = msg.sampleRate;
    if (typeof msg.frame_size === "number") FRAME_SIZE = msg.frame_size;
    if (typeof msg.overlapPercent === "number") OVERLAP_PERCENT = msg.overlapPercent;
    if (typeof msg.nfft === "number") NFFT = msg.nfft;

    if (typeof msg.chosenWindow === "string") CHOSEN_WINDOW = msg.chosenWindow;
    if (typeof msg.chosenMagnitude === "string") MAG_MODE = msg.chosenMagnitude;

    if (typeof msg.dbRef === "number") DB_REF = msg.dbRef;
    if (typeof msg.dbMin === "number") DB_MIN = msg.dbMin;
    if (typeof msg.dbMax === "number") DB_MAX = msg.dbMax;

    // reset FIFO on config changes (predictable)
    fifo = new Float32Array(0);
    return;
  }

  if (msg && msg.type === "paused") {
    PAUSED = !!msg.paused;
    return;
  }

  if (PAUSED) return;

  // Incoming audio stream chunk (Float32Array or transferable buffer)
  let chunk;
  if (msg instanceof Float32Array) chunk = msg;
  else chunk = new Float32Array(msg);

  fifo = appendFloat32(fifo, chunk);
  processFrames();
};
