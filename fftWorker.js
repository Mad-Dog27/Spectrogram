/* 
Creation: 04/07/2025
The purpose of this worker function is to handle the raw mic audio, it will resample


*/
let g =0;
let DEVICE_SAMPLE_RATE = 16000;
let SAMPLE_RATE = 16000;
let FRAME_SIZE = 128;
let CAPTURE_SIZE = 128;
let ratio =  SAMPLE_RATE/DEVICE_SAMPLE_RATE;
console.log("frf", ratio)
let lowerPower = 1;
let higherPower = 1;
let closestFrameSize = FRAME_SIZE;
let neededFrameSize = FRAME_SIZE / ratio;
let closestNeededFrameSize = neededFrameSize;
let NFFT = 4096
let expectedChunkTime = CAPTURE_SIZE/SAMPLE_RATE
let newAudioChunk = new Float32Array(FRAME_SIZE)

let currentBuffer = new Float32Array(FRAME_SIZE)
let prevBuffer = new Float32Array(FRAME_SIZE)

let overlapPercent = 0.25;
let overlap = Math.round(FRAME_SIZE * overlapPercent);
let PAUSED = false;

let chosenWindow = "blackman Harris"// rectangular, hamming, blackman Harris


// Called when main thread sends audio chunks
onmessage = function (e) {
    if (e.data.type == "config"){
        SAMPLE_RATE = e.data.sampleRate;
        DEVICE_SAMPLE_RATE = e.data.deviceSampleRate;
        FRAME_SIZE = e.data.frame_size
    } else if (e.data.type === "paused") {
        PAUSED = e.data.paused;
    }else {
        
        if (!PAUSED) {
           let start1 = performance.now();
        if (e.data.length*2 > NFFT) {NFFT = e.data.length*2}

        currentBuffer = new Float32Array(e.data);
        console.log(currentBuffer)
        const overlappedBuffer = addOverLap();
                console.log(overlappedBuffer)

        const preparedChunk = addZeroes(applyWindow(overlappedBuffer, FRAME_SIZE));
        const fftOutput = fft(preparedChunk)
        const halfLength = Math.floor(fftOutput.length / 2);
        //const halfFFT = fftOutput.slice(0, halfLength);
        //const flattened = new Float32Array(halfFFT.length * 2); // real + imag
            /*
        for (let i = 0; i < halfFFT.length; i++) {
            flattened[i * 2]     = halfFFT[i].real;
            flattened[i * 2 + 1] = halfFFT[i].imag;
        }*/
        const len = fftOutput.length;
        const magnitudes = new Float32Array(len);

        for (let i = 0; i < len; i++) {
        const real = fftOutput[i].real;
        const imag = fftOutput[i].imag;
        magnitudes[i] = Math.sqrt(real * real + imag * imag);
        }
      
        let start2 = performance.now();
        //console.log("fft time: ", start2 - start1)
        postMessage(magnitudes, [magnitudes.buffer]); // Pass along to next stage (fftWorker)

    }   
}
};

function addZeroes(frame) { // add zeroes to match nFFT value 
    let N = frame.length;
    if (N > NFFT) N = NFFT; // edge case, SHOULD NEVER HAPPEN

    if (N == NFFT) return frame; // no change needed
    const numZeroes = NFFT - N;
    const leftZeroes = Math.floor(numZeroes / 2); // zero padding to left of samples

    const paddedFrame = new Float32Array(NFFT);

    paddedFrame.set(frame, leftZeroes); // right zeropadding automatically done when creating new array ^^

    return paddedFrame;
}
function applyWindow(chunk, frameLength) {
    if (chosenWindow == "rectangular") { // no change
        return chunk;
    }

    if (chosenWindow == "hamming") {
        for (let n = 0; n < frameLength; n++) {
            chunk[n] = chunk[n] * (0.54 - 0.46 * Math.cos((2 * Math.PI * n) / (frameLength - 1)));
        }


    }
    if (chosenWindow == "blackman Harris") {
        for (let n = 0; n < frameLength; n++) {
            chunk[n] = chunk[n] * (0.35875 - 0.48829 * Math.cos((2 * Math.PI * n) / (frameLength - 1)) +
                0.14128 * Math.cos((4 * Math.PI * n) / (frameLength - 1)) -
                0.01168 * Math.cos((6 * Math.PI * n) / (frameLength - 1)));
        }

    }
    return chunk;

}
// addoverlap to a buffer using the previous buffer
function addOverLap() {
    const safeCurrent = new Float32Array(currentBuffer); // full copy
    const safePrev = new Float32Array(prevBuffer);       // full copy

  let newOverlap = Math.floor(overlapPercent * FRAME_SIZE);
  const newCurrentBuffer = new Float32Array(FRAME_SIZE);
  // If no previous buffer, just fill with zeros + current
  if (safePrev.length === 0) {
    newCurrentBuffer.set(safeCurrent.subarray(0, FRAME_SIZE - newOverlap), newOverlap);
    safePrev = newCurrentBuffer.slice(); // store for next round
    return newCurrentBuffer;
  }

  prevLength = safePrev.length
  currentLength = safeCurrent.length
  if (prevLength < currentLength) {
    console.log("PREV LE: ", prevLength)
    newOverlap = prevLength 
    console.log("PREV LE: ", newOverlap)

  }

  // Copy overlap from end of previous buffer to start of new buffer
  newCurrentBuffer.set(safePrev.slice(safePrev.length - newOverlap), 0);
  console.log(newCurrentBuffer)

  // Copy new audio into rest
  newCurrentBuffer.set(safeCurrent.slice(0, FRAME_SIZE - newOverlap), newOverlap);

  // Update prevBuffer for next call
/*
      console.log(safePrev.slice(safePrev.length - newOverlap))
  console.log(safeCurrent.slice(0, FRAME_SIZE - newOverlap))
  console.log(newCurrentBuffer)*/


    prevBuffer = newCurrentBuffer.slice();

  return newCurrentBuffer;
}


 function bitReverseIndex(index, bits) {
    let reversed = 0;
    for (let i = 0; i < bits; i++) {
        reversed <<= 1;
        reversed |= index & 1;
        index >>= 1;
    }
    return reversed;
}

function fft(input) {
    const N = input.length;
    const levels = Math.log2(N);

    if ((N & (N - 1)) !== 0) {
        throw new Error("Input length must be a power of 2");
    }

    // Create output array of complex numbers
    const output = new Array(N);
    for (let i = 0; i < N; i++) {
        const j = bitReverseIndex(i, levels);
        output[i] = {
            real: input[j],
            imag: 0
        };
    }

    for (let size = 2; size <= N; size <<= 1) {
        const halfSize = size / 2;
        const angleStep = (-2 * Math.PI) / size;

        for (let i = 0; i < N; i += size) {
            for (let j = 0; j < halfSize; j++) {
                const even = output[i + j];
                const odd = output[i + j + halfSize];

                const angle = angleStep * j;
                const twiddle = {
                    real: Math.cos(angle),
                    imag: Math.sin(angle)
                };

                const t = {
                    real: twiddle.real * odd.real - twiddle.imag * odd.imag,
                    imag: twiddle.real * odd.imag + twiddle.imag * odd.real
                };

                output[i + j] = {
                    real: even.real + t.real,
                    imag: even.imag + t.imag
                };

                output[i + j + halfSize] = {
                    real: even.real - t.real,
                    imag: even.imag - t.imag
                };
            }
        }
    }

    return output;
}