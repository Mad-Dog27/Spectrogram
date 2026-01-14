/* micWorker.js
   Purpose:
   - Receive device-sample-rate Float32 chunks from AudioWorklet/main
   - (Optional) low-pass filter when downsampling
   - Resample to target SAMPLE_RATE as a continuous stream (no dropped remainder)
   - Post resampled chunks to next stage (fftWorker)

   Notes:
   - Uses streaming polyline resampling (linear interp) with phase accumulator.
   - Includes optional FIR low-pass to reduce aliasing when downsampling.
*/

"use strict";

let DEVICE_SAMPLE_RATE = 48000;
let SAMPLE_RATE = 16000;

// Streaming input FIFO (device-rate)
let inBuf = new Float32Array(0);

// Resampler phase accumulator (in device-sample units)
let phase = 0;          // fractional index into inBuf
let step = 1;           // deviceSamplesPerTargetSample = DEVICE_SR / TARGET_SR

let PAUSED = false;

// FIR anti-alias (only needed if downsampling)
let USE_FIR = true;
let FIR_TAPS = 63;      // odd number
let firKernel = null;
let firHist = new Float32Array(0); // recent device samples for streaming FIR

function rebuildResampler() {
  step = DEVICE_SAMPLE_RATE / SAMPLE_RATE;
  // reset phase to avoid weirdness when changing rates
  phase = 0;

  if (USE_FIR && SAMPLE_RATE < DEVICE_SAMPLE_RATE) {
    // cutoff slightly below target Nyquist, expressed at device sample-rate
    const cutoffHz = 0.45 * SAMPLE_RATE; // conservative
    firKernel = makeLowpassFIR(cutoffHz, DEVICE_SAMPLE_RATE, FIR_TAPS);
    firHist = new Float32Array(FIR_TAPS - 1); // keep previous samples for convolution
  } else {
    firKernel = null;
    firHist = new Float32Array(0);
  }
}

function appendFloat32(a, b) {
  if (a.length === 0) return b;
  if (b.length === 0) return a;
  const out = new Float32Array(a.length + b.length);
  out.set(a, 0);
  out.set(b, a.length);
  return out;
}

// Windowed-sinc lowpass FIR (Hamming)
function makeLowpassFIR(cutoffHz, sampleRate, taps) {
  const h = new Float32Array(taps);
  const M = taps - 1;
  const fc = cutoffHz / sampleRate; // normalized (0..0.5)

  let sum = 0;
  for (let n = 0; n < taps; n++) {
    const k = n - M / 2;
    let sinc;
    if (k === 0) sinc = 2 * fc;
    else sinc = Math.sin(2 * Math.PI * fc * k) / (Math.PI * k);

    // Hamming window
    const w = 0.54 - 0.46 * Math.cos((2 * Math.PI * n) / M);
    h[n] = sinc * w;
    sum += h[n];
  }

  // normalize DC gain
  for (let n = 0; n < taps; n++) h[n] /= sum;
  return h;
}

// Streaming FIR: filter a device-rate block, maintaining history
function firFilterBlock(x) {
  // prepend history so convolution is valid at block start
  const xExt = appendFloat32(firHist, x);
  const y = new Float32Array(x.length);
  const taps = firKernel.length;
  const half = taps - 1; // history length

  for (let i = 0; i < x.length; i++) {
    let acc = 0;
    const center = i + half; // position in xExt aligned to output sample i
    for (let k = 0; k < taps; k++) {
      acc += xExt[center - k] * firKernel[k];
    }
    y[i] = acc;
  }

  // update history with last (taps-1) samples of xExt
  firHist = xExt.slice(xExt.length - (taps - 1));
  return y;
}

// Produce a chunk of resampled audio from inBuf.
// Returns Float32Array (may be empty if not enough samples).
function resampleFromBuffer(maxOut = 1024) {
  // Need at least 2 samples beyond phase for linear interpolation.
  // Also, if FIR enabled, inBuf has already been filtered blockwise.
  const out = new Float32Array(maxOut);
  let outCount = 0;

  while (outCount < maxOut) {
    const i0 = Math.floor(phase);
    const frac = phase - i0;
    const i1 = i0 + 1;

    if (i1 >= inBuf.length) break; // not enough input yet

    const s0 = inBuf[i0];
    const s1 = inBuf[i1];
    out[outCount++] = s0 + frac * (s1 - s0);

    phase += step;
  }

  // Drop consumed input samples to keep buffer small.
  // Keep one sample before i0 for continuity.
  const drop = Math.max(0, Math.floor(phase) - 1);
  if (drop > 0) {
    inBuf = inBuf.slice(drop);
    phase -= drop;
  }

  return outCount === maxOut ? out : out.slice(0, outCount);
}

rebuildResampler();

onmessage = (e) => {
  const msg = e.data;

  if (msg && msg.type === "config") {
    if (typeof msg.sampleRate === "number") SAMPLE_RATE = msg.sampleRate;
    if (typeof msg.deviceSampleRate === "number") DEVICE_SAMPLE_RATE = msg.deviceSampleRate;
    if (typeof msg.useFIR === "boolean") USE_FIR = msg.useFIR;
    if (typeof msg.firTaps === "number" && msg.firTaps >= 15 && (msg.firTaps % 2 === 1)) FIR_TAPS = msg.firTaps;

    // Reset buffers on config change (keeps behaviour predictable)
    inBuf = new Float32Array(0);
    phase = 0;
    rebuildResampler();
    return;
  }

  if (msg && msg.type === "paused") {
    PAUSED = !!msg.paused;
    return;
  }

  if (PAUSED) return;

  // Raw audio chunk (Float32Array or transferable ArrayBuffer)
  let chunk;
  if (msg instanceof Float32Array) chunk = msg;
  else chunk = new Float32Array(msg); // if we received a buffer

  // Optional FIR for anti-alias when downsampling
  if (firKernel) {
    chunk = firFilterBlock(chunk);
  }

  // Append to device-rate buffer
  inBuf = appendFloat32(inBuf, chunk);

  // Emit resampled stream chunks (small bursts)
  // You can tune maxOut: larger = fewer messages, smaller = lower latency
  while (true) {
    const out = resampleFromBuffer(1024);
    if (out.length === 0) break;

    // Transferable
    postMessage(out, [out.buffer]);
  }
};
