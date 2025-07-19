/* 
Creation: 04/07/2025
The purpose of this worker function is to handle the raw mic audio, it will resample


*/

let DEVICE_SAMPLE_RATE = 48000;
let SAMPLE_RATE = 16000;
let FRAME_SIZE = 128;
let PAUSED = false;
let CAPTURE_SIZE = 128;
let ratio =  SAMPLE_RATE/DEVICE_SAMPLE_RATE;

let closestFrameSize = FRAME_SIZE;
let neededFrameSize = FRAME_SIZE / ratio;
let closestNeededFrameSize = neededFrameSize;
let NFFT = 4096
let expectedChunkTime = CAPTURE_SIZE/16000

let newAudioChunk = new Float32Array(128)
let totalAudioChunk = null;
let chosenChunk = null;

let KERNEL = generateLowPassKernel(SAMPLE_RATE/2, DEVICE_SAMPLE_RATE, 101)
updateRequiredFrameSize();


// Called when main thread sends audio chunks
onmessage = function (e) {
    if (e.data.type === "config") {
    SAMPLE_RATE = e.data.sampleRate;
    DEVICE_SAMPLE_RATE = e.data.deviceSampleRate;
    FRAME_SIZE = e.data.frame_size
    closestFrameSize = FRAME_SIZE;
    ratio =  SAMPLE_RATE/DEVICE_SAMPLE_RATE;
    
    neededFrameSize = FRAME_SIZE / ratio;
    KERNEL = generateLowPassKernel(SAMPLE_RATE/2, DEVICE_SAMPLE_RATE, 101)
    updateRequiredFrameSize();
    


    } else if (e.data.type === "paused") {
        PAUSED = e.data.paused;
    } else { 
        if (!PAUSED) {
        
            let audioChunk = e.data

            if (newAudioChunk == null) {
               newAudioChunk = audioChunk
               return
            }
            newAudioChunk = appendBuffer(newAudioChunk, audioChunk);

            if (newAudioChunk.length >= neededFrameSize) {
                totalAudioChunk = newAudioChunk.slice(0, neededFrameSize)
                newAudioChunk = null;
                if (SAMPLE_RATE != DEVICE_SAMPLE_RATE) {
                let filteredChunk = applyFIRFilter(totalAudioChunk, KERNEL)
                //let resampledAudioChunk = downsample(filteredChunk, 3)
                let resampledAudioChunk = resampleMicBuffer(filteredChunk);
                    chosenChunk = resampledAudioChunk;
                } else {
                    chosenChunk = totalAudioChunk;
                }
                let thunk = totalAudioChunk
                postMessage(chosenChunk, [chosenChunk.buffer]); // Pass along to next stage (fftWorker)
                   
                totalAudioChunk = null;

             } else {return;} 
      
    
    let start2 = performance.now();
    //postMessage({ type: "print", currentBuffer});


    //console.log("Sampling time: ", start2 - start1)
    }
}
};

function appendBuffer(buffer1, buffer2) {
  const tmp = new Float32Array(buffer1.length + buffer2.length);
  tmp.set(buffer1, 0);
  tmp.set(buffer2, buffer1.length);
  return tmp;
}

function updateRequiredFrameSize() {
            let lowerPower = 1;
            let higherPower = 1;    
            while (lowerPower * 2 < neededFrameSize) { lowerPower <<= 1; }
            higherPower = lowerPower * 2;
            lowerPower >>= 1;
            closestNeededFrameSize = higherPower;
            higherPower = 1;
            lowerPower = 1;
            // calculating closest frame size to the power of 2 to the users desired. This is redundant right now
            while (lowerPower < FRAME_SIZE) { lowerPower <<= 1; }
            higherPower = lowerPower;
            lowerPower >>= 1;
            closestFrameSize = higherPower;
            //analyser.fftSize = closestNeededFrameSize; 
            
            if (closestFrameSize > NFFT) { 
                NFFT = closestFrameSize * 2; //frequency domain amount zeroes and values aquired through fft
            }
            
}
function resampleMicBuffer(buffer) {//This works
    //Buffer is going to be ratio times larger then need be, this is becuase of interpolation, lower freqs skip a certain
    // amount of samples, but still needs to be same size, therefore original bufffer needs to be larger
    //console.log(buffer)
    //const ratio = DEVICE_SAMPLE_RATE / SAMPLE_RATE; // deviceFS / sampleFREQ (user chosen)
    let newBuffer = new Float32Array(FRAME_SIZE);
    if (DEVICE_SAMPLE_RATE == SAMPLE_RATE) { return buffer; } // if no resampling needed
    for (let i = 0; i < FRAME_SIZE; i++) { //newLength should be buffers original size(desired size)
        // by grabbing every ratio sample, means sample freq is lowered by ratio
        let index = i / ratio;
        let lowerIndex = Math.floor(index);
        let upperIndex = Math.ceil(index);
        let fraction = index - lowerIndex; // if not a clean ratio, ie decimal
        if (upperIndex < buffer.length) { //Interpulation (if needed)
            newBuffer[i] = buffer[lowerIndex] * (1 - fraction) + buffer[upperIndex] * fraction; // Linear interpolation
        } else {
            newBuffer[i] = buffer[lowerIndex]; // Edge case handling
        }
    }
    //console.log(newBuffer[newBuffer.length-1])
    return newBuffer;
}


function generateLowPassKernel(cutoffFreq, sampleRate, kernelSize) {
    const kernel = new Float32Array(kernelSize);
    const fc = cutoffFreq / sampleRate; // normalized cutoff (0 < fc < 0.5)
    const M = kernelSize - 1;
    const PI = Math.PI;

    for (let n = 0; n < kernelSize; n++) {
        if (n === M / 2) {
            kernel[n] = 2 * fc;
        } else {
            const x = n - M / 2;
            kernel[n] = Math.sin(2 * PI * fc * x) / (PI * x);
        }

        // Apply Hamming window
        kernel[n] *= 0.54 - 0.46 * Math.cos((2 * PI * n) / M);
    }

    return kernel;
}

function applyFIRFilter(input, kernel) {
    const output = new Float32Array(input.length);
    const half = Math.floor(kernel.length / 2);

    for (let i = 0; i < input.length; i++) {
        let acc = 0;
        for (let j = 0; j < kernel.length; j++) {
            const idx = i - j + half;
            if (idx >= 0 && idx < input.length) {
                acc += input[idx] * kernel[j];
            }
        }
        output[i] = acc;
    }

    return output;
}

function downsample(input, factor) {
    const outputLength = Math.floor(input.length / factor);
    const output = new Float32Array(outputLength);

    for (let i = 0; i < outputLength; i++) {
        output[i] = input[i * factor];
    }

    return output;
}