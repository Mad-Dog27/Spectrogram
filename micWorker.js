/* 
Creation: 04/07/2025
The purpose of this worker function is to handle the raw mic audio, it will resample


*/

let DEVICE_SAMPLE_RATE = 16000;
let SAMPLE_RATE = 16000;
let FRAME_SIZE = 128;
let PAUSED = false;
let CAPTURE_SIZE = 128;
let ratio =  SAMPLE_RATE/DEVICE_SAMPLE_RATE;
let lowerPower = 1;
let higherPower = 1;
let closestFrameSize = FRAME_SIZE;
let neededFrameSize = FRAME_SIZE / ratio;
let closestNeededFrameSize = neededFrameSize;
let NFFT = 4096
let expectedChunkTime = CAPTURE_SIZE/SAMPLE_RATE
let newAudioChunk = new Float32Array(128)


updateRequiredFrameSize();


// Called when main thread sends audio chunks
onmessage = function (e) {
    if (e.data.type === "config") {
    SAMPLE_RATE = e.data.sampleRate;
    DEVICE_SAMPLE_RATE = e.data.deviceSampleRate;
    updateRequiredFrameSize();

    } else if (e.data.type === "pause") {
        PAUSED = e.data.paused;
    } else { 
    if (!PAUSED) {
    let start1 = performance.now();

    let i = 0;   
    let prevTime = 0;
    let startTime = performance.now();
    let timePassed=0;
    while (i < ratio) {
        if (timePassed > expectedChunkTime) {
            const audioChunk = new Float32Array(e.data); 
            newAudioChunk = appendBuffer(audioChunk, newAudioChunk)
            i++;
            timePassed = 0;
            prevTime = 0;
        }
        thisIterationTime = performance.now() - startTime - prevTime; 
        prevTime = thisIterationTime;
        timePassed += thisIterationTime;
    }
    resampledAudioChunk = resampleMicBuffer(newAudioChunk);
 
  newAudioChunk.set(0)
    let start2 = performance.now();
    //console.log("Sampling time: ", start2 - start1)
    postMessage(resampledAudioChunk, [resampledAudioChunk.buffer]); // Pass along to next stage (fftWorker)
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

    //const ratio = DEVICE_SAMPLE_RATE / SAMPLE_RATE; // deviceFS / sampleFREQ (user chosen)
    let newBuffer = new Float32Array(FRAME_SIZE);
    if (DEVICE_SAMPLE_RATE == SAMPLE_RATE) { return buffer; } // if no resampling needed
    for (let i = 0; i < FRAME_SIZE; i++) { //newLength should be buffers original size(desired size)
        // by grabbing every ratio sample, means sample freq is lowered by ratio
        let index = i * ratio;
        let lowerIndex = Math.floor(index);
        let upperIndex = Math.ceil(index);
        let fraction = index - lowerIndex; // if not a clean ratio, ie decimal

        if (upperIndex < buffer.length) { //Interpulation (if needed)
            newBuffer[i] = buffer[lowerIndex] * (1 - fraction) + buffer[upperIndex] * fraction; // Linear interpolation
        } else {
            newBuffer[i] = buffer[lowerIndex]; // Edge case handling
        }
    }
    
    return newBuffer;
}