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
let lowerPower = 1;
let higherPower = 1;
let closestFrameSize = FRAME_SIZE;
let neededFrameSize = FRAME_SIZE / ratio;
let closestNeededFrameSize = neededFrameSize;
let NFFT = 8192
let expectedChunkTime = CAPTURE_SIZE/SAMPLE_RATE
let newAudioChunk = new Float32Array(FRAME_SIZE)

let currentBuffer = new Float32Array(FRAME_SIZE)
let prevBuffer = new Float32Array(FRAME_SIZE)

let overlapPercent = 0;
let overlap = Math.round(FRAME_SIZE * overlapPercent);
let PAUSED = false;

let CHOSEN_MAGNITUDE_SCALE = "magnitude"

let OVERLAP_PERCENT = 0.25
let CHOSEN_WINDOW = "blackman Harris"

let movingAvg = new Float32Array(128)

// Called when main thread sends audio chunks
onmessage = function (e) {
    if (e.data.type == "config"){
        if ("sampleRate" in e.data) SAMPLE_RATE = e.data.sampleRate;
        if ("deviceSampleRate" in e.data) DEVICE_SAMPLE_RATE = e.data.deviceSampleRate;
        if ("frame_size" in e.data) FRAME_SIZE = e.data.frame_size;
        if ("overlapPercent" in e.data) OVERLAP_PERCENT = e.data.overlapPercent;
        if ("chosenWindow" in e.data) CHOSEN_WINDOW = e.data.chosenWindow;
        if ("chosenMagnitude" in e.data) CHOSEN_MAGNITUDE_SCALE = e.data.chosenMagnitude
        console.log("fs (target):", SAMPLE_RATE, ", fs (device):",DEVICE_SAMPLE_RATE, ", frame size:", FRAME_SIZE,", overlap:", OVERLAP_PERCENT*100,"%, Window:", CHOSEN_WINDOW, ", Magnitude:",CHOSEN_MAGNITUDE_SCALE)

    } else if (e.data.type === "paused") {
        PAUSED = e.data.paused;
    }else {
        
        if (!PAUSED) {
            let chosenValues;

            currentBuffer = new Float32Array(e.data);
            const overlappedBuffer = addOverLap(currentBuffer)
            const chunk = addZeroes(applyWindow(overlappedBuffer, FRAME_SIZE));//Applying a window AND zero padding, the function above defaults to rectangular window

            const thisChunk = fft(chunk)

            if (CHOSEN_MAGNITUDE_SCALE == "magnitude") {
                const data = computeMagnitudeAndDB(thisChunk, 1, false);
                let dataMagnitude = data.map(bin => bin.magnitude);
                chosenValues = dataMagnitude.slice(0, NFFT / 2);
                
                /*currentValues2 = currentBuffer.map((val, i) => 0.2*val * movingAvg[i])
                
                chosenValues = currentValues.map((val, i) => val*0.8) + currentValues2
                movingAvg = new Float32Array(currentValues);
                */
              

            } else {
                const data = computeMagnitudeAndDB(thisChunk, 1, true);
                let datadB = data.map(bin => bin.dB);
                chosenValues = datadB.slice(0, NFFT / 2)
            }
            let fftOutput = new Float32Array(chosenValues)

            postMessage({ type: "print", currentBuffer});
            postMessage(fftOutput, [fftOutput.buffer]);
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
function applyWindow(unwindowChunk, frameLength) {
    let chunk = new Float32Array(unwindowChunk)
    if (CHOSEN_WINDOW == "rectangular") { // no change
        return chunk;
    }

    if (CHOSEN_WINDOW == "hamming") {
        for (let n = 0; n < frameLength; n++) {
            chunk[n] = unwindowChunk[n] * (0.54 - 0.46 * Math.cos((2 * Math.PI * n) / (frameLength - 1)));
        }


    }
    if (CHOSEN_WINDOW == "blackman Harris") {
        for (let n = 0; n < frameLength; n++) {
            chunk[n] = unwindowChunk[n] * (0.35875 - 0.48829 * Math.cos((2 * Math.PI * n) / (frameLength - 1)) +
                0.14128 * Math.cos((4 * Math.PI * n) / (frameLength - 1)) -
                0.01168 * Math.cos((6 * Math.PI * n) / (frameLength - 1)));
        }

    }
    return chunk;

}
// addoverlap to a buffer using the previous buffer
function addOverLap() {
    if (OVERLAP_PERCENT > 0) {
    const safeCurrent = new Float32Array(currentBuffer); // full copy
    const safePrev = new Float32Array(prevBuffer);       // full copy

  let newOverlap = Math.floor(OVERLAP_PERCENT * FRAME_SIZE);
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
    newOverlap = prevLength 
  }

  newCurrentBuffer.set(safePrev.slice(safePrev.length - newOverlap), 0);
  newCurrentBuffer.set(safeCurrent.slice(0, FRAME_SIZE - newOverlap), newOverlap);

 


    prevBuffer = newCurrentBuffer.slice();

  return newCurrentBuffer;
    }
    return currentBuffer
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
function computeMagnitudeAndDB(fftResult, REF, useDB) {
    return fftResult.map(({ real, imag }) => {
        const magnitude = Math.sqrt(real ** 2 + imag ** 2);
        let dB = 0;
        if ((magnitude > 0) && (useDB)) {
        dB = useDB ? 20 * Math.log10(magnitude / REF) : -100;
        } 
        return { real, imag, magnitude, dB };
    });
} 