/* 
Creation: 04/07/2025
The purpose of this worker function is to handle the raw mic audio, it will resample


*/

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

let chosenWindow = "blackman Harris"// rectangular, hamming, blackman Harris


// Called when main thread sends audio chunks
onmessage = function (e) {
    if (e.data.type == "config"){

    } else {
        currentBuffer = new Float32Array(e.data);
        const overlappedBuffer = addOverLap();
        const preparedChunk = addZeroes(applyWindow(overlappedBuffer, FRAME_SIZE));
        const fftOutput = fft(preparedChunk)
        
        const flattened = new Float32Array(fftOutput.length * 2); // real + imag

        for (let i = 0; i < fftOutput.length; i++) {
            flattened[i * 2]     = fftOutput[i].real;
            flattened[i * 2 + 1] = fftOutput[i].imag;
        }
        postMessage(flattened, [flattened.buffer]); // Pass along to next stage (fftWorker)

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
function addOverLap() { //WRONG OVERLAP SHOUDLNT INCREASE FRAMESIZE
    
    const prevLength = prevBuffer.length;
    //new overlap in terms of samples
    const newOverlap = Math.floor(overlapPercent * FRAME_SIZE);
    
    let newCurrentBuffer = new Float32Array(FRAME_SIZE)
    if (prevLength == 0) { // if no previous buffer, then no overlap applied 
        newCurrentBuffer.set(currentBuffer.subarray(0, newCurrentBuffer.length - newOverlap), newOverlap);
        return newCurrentBuffer;
    }

    for (let i = 0; i < newOverlap; i++) { // add end overlap amount of samples of previous buffer to the start of the new buffer
        newCurrentBuffer[i] = prevBuffer[prevLength - newOverlap + i]
    }
    // add the rest of the buffer
    newCurrentBuffer.set(currentBuffer.subarray(0, newCurrentBuffer.length - newOverlap), newOverlap)

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