/* READ ME
Purpose of the file is to produce an interactable spectrogram for users to view how their pronounications look,
while providing the ability to alter the visuals of it. 
Current uses: -Audio File Input
              -Microphone Input
              -Recording and Save mircophone Input
              -Mel Scaling
              -DB magnitude - Not working
              -Colour Scheme
              -Window Type

              OUTPUTS:
              -Spectrogram Output
              -Waveform Output
              -Audio Output - Recorded Samples sound weird (sampling frequency problem)
            
To do:

Make website look way nicer
Fix Colour Scheme
Fix standered units output from FFT
Fix DB
Fix Mel
Replace all toggles with good toggles eg if else

NOTE: may crash your browser :) - Blame timeGraph() function
*/

const audioContext = new (window.AudioContext || window.webkitAudioContext)();

//Button INputs for inputing Data
const audioFileInput = document.getElementById('audioFileInput');
const processAgainInput = document.getElementById('processAgain');
const micButtonInput = document.getElementById('Mic');
const recordValue = document.getElementById('recordValue')
const playRecordButton = document.getElementById('playRecord')
//Button Inputs for altering graphs
const melButtonInput = document.getElementById("melToggle")
const isMel = document.getElementById('isMel')
const recordButtonInput = document.getElementById('recordToggle')
const windowSelect = document.getElementById('windowSelect');
const colourSchemeSelect = document.getElementById('colourSelect');
const magnitudeSelect = document.getElementById('magnitudeSelect')

const widthSlider = document.getElementById('widthSlider');
const widthSliderValue = document.getElementById('widthSliderValue');

const sampleFreqSlider = document.getElementById('sampleFreqSlider');
const sampleFreqSliderValue = document.getElementById('sampleFreqSliderValue');

const overlapPercSlider = document.getElementById('overlapPercSlider');
const overlapPercSliderValue = document.getElementById('overlapPercSliderValue');
//IDs for all three canvas's: Canvas2D, CanvasSpectrum, timeCanvas
const canvas2D = document.getElementById("canvas2D");
const ctx2D = canvas2D.getContext("2d");
const canvasSpectrum = document.getElementById("canvasSpectrum");

const canvasAxis = document.getElementById('canvasAxis');
const ctxAxis = canvasAxis.getContext('2d');

const ctxSpectrum = canvasSpectrum.getContext("2d");
const timeCanvas = document.getElementById('timeCanvas')
const ctxTime = timeCanvas.getContext("2d")


//Global Constants
const FRAMESIZE = 1024; //time domain amount of samples taken
const nFFT = 2048; //frequency domain amount zeroes and values aquired through fft
let overlapPercent = 0.25;
let overlap = FRAMESIZE * overlapPercent;
const SPEED = 1;
let SAMPLEFREQ = 16000;
//Global Variables
let SCALE = 3;
let SENS = 1;
let CONTRAST = 0;
let recordOn = false;
let storedBuffer = [];
let audioBuffer;
let fileUpload = null;
let filePlaying = null;
let finshedRT = null;
let micOn = null;
let mel = null;
let melOn = false;
let melProcessed = false;

let REFERENCE = 2 ^ 15;
let WIDTH = 0.7;    //0.7
let HEIGHT = 0.7;  //0.49  DOESNT EFFECT - No point
let RecordProcessing = false;
chosenWindow = "blackman Harris"// rectangular, hamming, blackman Harris
chosenColourScheme = 'greyScale'
chosenMagnitudeScale = "magnitude"

canvasSpectrum.width = window.innerWidth * WIDTH - 2;  // 70% of screen width minus borders
canvasSpectrum.height = window.innerHeight * HEIGHT - 2;
canvasAxis.width = canvasSpectrum.width + 40;
canvasAxis.height = canvasSpectrum.height;

timeCanvas.width = window.innerWidth * (WIDTH * 3) - 2;  // 70% of screen width minus borders
timeCanvas.height = window.innerHeight * 0.45 - 2;  // 45% of screen height minus borders
drawAxisLabel();



//Adds an event listnener for the audioFileInput button, when the input file is changed the function will run after the file is changed.
audioFileInput.addEventListener('change', async (event) => { //AUDIO FILE INPUT
    const file = event.target.files[0]; // file list of all the inputed files, [0] means only one file will be used
    if (!file) return;
    if (micOn) { console.log("Microphone is still recording"); return; }
    fileUpload = true;
    const arrayBuffer = await file.arrayBuffer(); //Read the file as an ArrayBuffer which is a binary representation of the audio file to use in the next line
    audioBuffer = await audioContext.decodeAudioData(arrayBuffer);//Uses the binary version to create an audio buffer

    //mel = melScale()
    // Process the audio buffer (e.g., generate a spectrogram)
    filePlaying = true;
    processAudioBuffer(audioBuffer);
});


processAgainInput.addEventListener('click', () => {//  PROCCESS AGAIN
    //Event listener to process the audio file again, will happen on click
    if (!fileUpload) { console.log("No input file selected"); return; }
    if (filePlaying) { console.log("File is currently playing"); return; }
    if (micOn) { console.log("Mircophone Audio is currently playing"); return; }
    filePlaying = true;

    console.log("proccessing again")
    processAudioBuffer(audioBuffer);

})
micButtonInput.addEventListener('click', () => {
    //Event listener for Microphone input button, occurs on click, button acts as a toggle 
    if (filePlaying) { console.log("File Audio is currently playing"); return; }
    if (micOn) { // If mic already on, toggle off, NOTE: There is a better way to toggle, EDIT LATER
        micOn = false;
        audioContext.suspend();//Pause the audiocontext from capturing data
    } else { // If mic not on alreay, toggle on
        console.log("Listening to Mic")

        micOn = true;
        audioContext.resume()//Allow audio context to capture data
        getMicData();
    }
})
melButtonInput.addEventListener('click', () => {

    if (melOn) {
        melOn = null;
        drawAxisLabel()
        isMel.textContent = "off";
    } else {
        melOn = true;
        mel = melScale();
        drawAxisLabel()
        isMel.textContent = "on";
    }
    console.log("Mel toggled ", melOn)

})

magnitudeSelect.addEventListener('change', (event) => {
    //Changing the window value, between Rectangular, hamming and blackman-harris
    chosenMagnitudeScale = event.target.value;
    console.log(`Window function selected: ${chosenMagnitudeScale}`);
});
windowSelect.addEventListener('change', (event) => {
    //Changing the window value, between Rectangular, hamming and blackman-harris
    chosenWindow = event.target.value;
    console.log(`Window function selected: ${chosenWindow}`);
});
colourSchemeSelect.addEventListener('change', (event) => {
    //Changing the colour scheme, between grey scale and heated metal - WORK IN PROGRESS
    chosenColourScheme = event.target.value;
    console.log(`Colour Scheme Selected: ${chosenColourScheme}`);
});

widthSlider.addEventListener('input', () => {
    //Sensitivity slider display and use, on new input
    WIDTH = widthSlider.value; //Storing new value in SENS
    canvasSpectrum.width = window.innerWidth * WIDTH - 2;  // 70% of screen width minus borders
    timeCanvas.width = window.innerWidth * (WIDTH * 3) - 2;  // 70% of screen width minus borders

    widthSliderValue.textContent = WIDTH; // Update the display
    console.log(`Width: ${WIDTH}`);
});

sampleFreqSlider.addEventListener('input', () => {//Function to update the INputed Sampling freq, this will improved freqeuency resolution but only untill you reach the original inputed frequency
    SAMPLEFREQ = sampleFreqSlider.value; //
    sampleFreqSliderValue.textContent = SAMPLEFREQ; // Update the display
    mel = melScale();
    drawAxisLabel();

    console.log(`FS: ${SAMPLEFREQ}`);
});

overlapPercSlider.addEventListener('input', () => {//Function to update the INputed Sampling freq, this will improved freqeuency resolution but only untill you reach the original inputed frequency
    overlapPercent = overlapPercSlider.value; //
    overlapPercSliderValue.textContent = overlapPercent; // Update the display
    overlap = FRAMESIZE * overlapPercent;

});

recordButtonInput.addEventListener('click', () => {
    //Event listener to process the audio file again, will happen on click

    if (recordOn) {
        recordOn = false;
        recordValue.textContent = "no"
        //processRecording()
    } else {
        recordOn = true;
        recordValue.textContent = "yes"
    }

})
playRecordButton.addEventListener('click', () => {
    //Event listener to process the audio file again, will happen on click
    processRecording();

})





function processRecording() {
    audioBuffer = null;
    let combinedBuffer = combineBuffers(storedBuffer);

    audioBuffer = audioContext.createBuffer(
        1, // Number of channels
        combinedBuffer.length, // Length of the buffer
        SAMPLEFREQ // Sample rate
    );

    // Copy the combined buffer into the AudioBuffer
    audioBuffer.getChannelData(0).set(combinedBuffer);

    // Create a source and play the buffer
    const source = audioContext.createBufferSource();
    source.buffer = audioBuffer;
    source.connect(audioContext.destination);
    source.start(0);
    RecordProcessing = true;

    processAudioBuffer(audioBuffer);

}

async function getMicData() {
    try {
        if (mel == null) { }//mel = melScale(); }
        // Requesting microphone acces and waiting, then creating a mediastream from it
        const stream = await navigator.mediaDevices.getUserMedia({ audio: true });
        const mediaStreamSource = audioContext.createMediaStreamSource(stream);
        const deviceSampleRate = audioContext.sampleRate;

        // Create an AnalyserNode for real-time frequency domain analysis
        const analyser = audioContext.createAnalyser();
        analyser.fftSize = nFFT;
        mediaStreamSource.connect(analyser);

        // Array to store time-domain or frequency-domain data
        const dataArray = new Float32Array(analyser.frequencyBinCount);
        let chunkIndex = 0;
        // Function to process and visualize the microphone input 
        function processMicInput() {
            //New array/Buffer to store the audio samples
            let timeDomainBuffer = new Float32Array(FRAMESIZE);
            //analyser.sampleRate = 16000;
            //console.log(analyser.sampleRate)

            /*In the audio file processing, I store the audio data myself, but in this case I am not able to access
            the audio buffer, so am making use of getFloatTimeDomainData which stores the audio data from the analyser
            int my timeDomainBuffer. 
            NOTE: getFloatFrequencyData could also be used to obtain the frequency magnitudes But I prefer to use my maths*/
            analyser.getFloatTimeDomainData(timeDomainBuffer);

            const resampledTimeDomainBuffer = resampleBuffer(
                timeDomainBuffer,
                deviceSampleRate,
                SAMPLEFREQ
            )
            const chunk = addZeroes(applyWindow(resampledTimeDomainBuffer));//Applying a window AND zero padding, the function above defaults to rectangular window
            const result = fft(chunk); //FAST FOURIER TRANSFORM

            if (chosenMagnitudeScale == "magnitude") {
                const dataMagnitude = result.map(bin => bin.magnitude);
                chosenValues = dataMagnitude.slice(0, nFFT / 2)
            } else {
                const datadB = result.map(bin => bin.dB);
                chosenValues = datadB.slice(0, nFFT / 2)
            }
            //drawVisual(analyser)
            createSpectrum(chosenValues);
            createMovingSpectrogram(chosenValues);

            timeGraph(timeDomainBuffer)
            if (recordOn) {
                if (!storedBuffer[chunkIndex]) {
                    storedBuffer[chunkIndex] = resampledTimeDomainBuffer;
                } else {
                    console.warn("Chunk index already filled, potential overwrite detected.");
                }
                //storedBuffer[chunkIndex] = timeDomainBuffer;
                //storedBuffer[chunkIndex] = normalizeChunk(timeDomainBuffer, 0.9);
                chunkIndex++;
            }
            if (micOn) {
                requestAnimationFrame(processMicInput); //Keep calling proccessMicINput
            } else {
                mediaStreamSource.disconnect();
                analyser.disconnect();
                stream.getTracks().forEach(track => track.stop()); // Stops the media stream tracks
                console.log('Stream processing stopped.');
            }

        }

        // Start processing
        processMicInput();
    } catch (error) {
        console.error("Error accessing microphone:", error); //Might be better then using console.log
    }
}

function combineBuffers(buffers) {
    // Calculate total length, subtracting overlap for transitions
    let size = 0;
    buffers.forEach((chunk, index) => {
        size += chunk.length;
    })

    const totalLength = size;
    const combinedBuffer = new Float32Array(totalLength);

    let offset = 0; // Starting position in the combinedBuffer
    buffers.forEach((buffer, index) => {
        for (let i = 0; i < buffer.length; i++) {
            // Add the sample to the combinedBuffer, taking overlap into account
            combinedBuffer[offset + i] += buffer[i];
        }
        // Move the offset forward
        offset += buffer.length;
    });
    return combinedBuffer;
}
async function resampleAudio(audioBuffer) {
    const originalSampleFreq = audioBuffer.sampleRate;

    if (originalSampleFreq == SAMPLEFREQ) {
        return audioBuffer;
    }
    console.log("NOTRETURN")
    const offlineContext = new OfflineAudioContext(
        audioBuffer.numberOfChannels,
        Math.ceil(audioBuffer.length * (SAMPLEFREQ / originalSampleFreq)),
        SAMPLEFREQ
    );
    const source = offlineContext.createBufferSource();
    source.buffer = audioBuffer;
    source.connect(offlineContext.destination);

    source.start(0);

    const newBuffer = await offlineContext.startRendering();
    return newBuffer;
}

function resampleBuffer(buffer, deviceRate, targetRate) {
    const resampleRatio = deviceRate / targetRate;
    const newLength = Math.floor(buffer.length / resampleRatio);
    const resampledBuffer = new Float32Array(newLength);

    for (let i = 0; i < newLength; i++) {
        const index = i * resampleRatio;
        const lower = Math.floor(index);
        const upper = Math.ceil(index);

        if (upper < buffer.length) {
            const weight = index - lower;
            resampledBuffer[i] =
                (1 - weight) * buffer[lower] + weight * buffer[upper];
        } else {
            resampledBuffer[i] = buffer[lower];
        }
    }

    return resampledBuffer;
}


async function processAudioBuffer(audioBuffer) {
    //Outputting the Audiobuffer characterestics
    console.log('Audio Buffer:', audioBuffer);
    console.log('Sample Rate:', audioBuffer.sampleRate);
    console.log('Number of Channels:', audioBuffer.numberOfChannels);
    console.log('Duration (s):', audioBuffer.duration);
    console.log(audioBuffer.length);
    audioBuffer = await resampleAudio(audioBuffer);
    //Exectuting the FFT for file input
    const source = audioContext.createBufferSource();
    source.buffer = audioBuffer;
    source.connect(audioContext.destination);

    const analyser = audioContext.createAnalyser(); // Analyser node creation, helps extract information like frequencies or amplitudes
    source.connect(analyser); // Connect the audio source to the analyser
    analyser.connect(audioContext.destination); // Connects the analyser to the speakers
    executeFFTWithSync(audioBuffer, source, analyser);
}


function executeFFTWithSync(audioBuffer, source, analyser) {
    /* This fucntion deals with handling the processing, mathmatical operations and display of the audioBuffer,
    
    NOTE: Origninally it was only meant to deal with the FFT Maths but due to the nature of audio process, 
    much more needed to be added. In the future I would like to clean this function up 16/01/2025
    */
    //Dealing with the new AudioContext and audioBuffer and connecting the source to it

    //const audioContext = new AudioContext();


    /*For the next step, the audioBuffer needs to be cut into respective chunks for audioProccessing, 
    first the chunk characteristics need to be determined*/
    const sampleRate = audioBuffer.sampleRate;
    const effectiveChunkSize = FRAMESIZE - overlap; //Subracting overlap, as that portion of chunk is accounted for in the next
    const chunkDuration = effectiveChunkSize / sampleRate; // Duration of one chunk in seconds
    // Slice the audio into chunks
    const chunks = sliceIntoChunks(audioBuffer); // NOTE: sliceIntoChunks function also applies window
    const numChunks = chunks.length;
    let chosenValues;
    let result = [];
    let currentChunkIndex = 0;
    //mel = melScale();
    audioContext.resume()
    // Start audio playback
    source.start(0);
    const startTime = audioContext.currentTime;

    // Function to update the spectrogram
    function updateSpectrogram() {
        // Calculate the expected chunk index based on current playback time
        const elapsedTime = audioContext.currentTime - startTime;
        const expectedChunkIndex = Math.floor(elapsedTime / chunkDuration);
        // Process and render all chunks up to the expected index
        while (currentChunkIndex <= expectedChunkIndex && currentChunkIndex < numChunks) { //If chunk processing is slower then expected and chunk index has not exceded total number of chunks
            const chunk = addZeroes(chunks[currentChunkIndex]); //Zero Padding

            timeGraph(chunks[currentChunkIndex]);
            // Fast Fourier transform of current chunk, cutting results in half then converting to magnitude

            result = fft(chunk);
            //result[currentChunkIndex] = result[currentChunkIndex].slice(0, FRAMESIZE / 2);

            //Updating all Graphs 
            if (chosenMagnitudeScale == "magnitude") {
                const dataMagnitude = result.map(bin => bin.magnitude);
                chosenValues = dataMagnitude.slice(0, nFFT / 2)
            } else {
                const datadB = result.map(bin => bin.dB);
                chosenValues = datadB.slice(0, nFFT / 2)


            }
            //drawVisual(analyser)
            createSpectrum(chosenValues);
            createMovingSpectrogram(chosenValues);

            currentChunkIndex++; //incrementing chunk index
        }

        // Continue updating the spectrogram if there are more chunks to process
        if (currentChunkIndex < numChunks) {
            requestAnimationFrame(updateSpectrogram);
        } else {
            filePlaying = null;
        }
    }

    // Start updating the spectrogram
    requestAnimationFrame(updateSpectrogram);
}

function fft(input) {
    const N = input.length; //Assume N is of size of power 2, ie (2^n = N)
    if (N <= 1) return [{ real: input[0], imag: 0, magnitude: 0, dB: -Infinity }];

    if ((N & (N - 1)) !== 0) {
        console.log("Input array length must be a power of 2.");
        return [];
    }
    //Split input into evens and odds then pass them back into function, untill size is 1
    const even = fft(input.filter((_, i) => i % 2 === 0));
    const odd = fft(input.filter((_, i) => i % 2 !== 0));


    const magnitude = Array(N)
    const combined = Array(N).fill(0).map(() => ({ real: 0, imag: 0, magnitude: 0, dB: -Infinity }));

    for (let k = 0; k < N / 2; k++) {
        const angle = (-2 * Math.PI * k) / N;
        const twiddle = complexExp(angle);
        /*const twiddle = {
            real: Math.cos(angle),
            imag: Math.sin(angle)
        }*/
        const t = {
            real: twiddle.real * odd[k].real - twiddle.imag * odd[k].imag,
            imag: twiddle.real * odd[k].imag + twiddle.imag * odd[k].real
        }

        combined[k] = {
            real: even[k].real + t.real,
            imag: even[k].imag + t.imag
        }
        combined[k + N / 2] = {
            real: even[k].real - t.real,
            imag: even[k].imag - t.imag
        }

        combined[k].magnitude = Math.sqrt(combined[k].real ** 2 + combined[k].imag ** 2);
        combined[k].dB = 20 * Math.log10(combined[k].magnitude);

        combined[k + N / 2].magnitude = Math.sqrt(
            combined[k + N / 2].real ** 2 + combined[k + N / 2].imag ** 2
        );
        combined[k + N / 2].dB = 20 * Math.log10(combined[k + N / 2].magnitude / 1);
        //magnitude[k] = Math.sqrt(combined[k].real ** 2 + combined[k].imag ** 2);
    }
    //console.log(Math.max(...combined.map(c => c.magnitude)));
    return combined;
}
function sliceIntoChunks(audioBuffer) {
    let tempOverlap = overlap;

    if (RecordProcessing) {//DOUBLE CHECK THIS
        overlap = 0;
    }
    const monoChannel = audioBuffer.getChannelData(0);
    if (audioBuffer.numberOfChannels >= 2) {
        const samples1 = audioBuffer.getChannelData(1);
        for (let j = 0; j < audioBuffer.length; j++) {
            monoChannel[j] = (monoChannel[j] + samples1[j]) / 2
        }
    }
    const numChunks = Math.floor((audioBuffer.length - FRAMESIZE) / (FRAMESIZE - overlap)) + 1;
    const chunks = [];


    for (let i = 0; i < numChunks; i++) {
        const start = i * (FRAMESIZE - overlap);
        const end = start + FRAMESIZE;
        const chunk = monoChannel.slice(start, end);
        const windowedChunk = applyWindow(chunk);
        //const chunk = audioBuffer.getChannelData(0).slice(offset, offset + bufferLength);
        chunks.push(windowedChunk);
    }
    overlap = tempOverlap;
    return chunks
}

function applyWindow(chunk) {
    if (chosenWindow == "rectangular") {
        return chunk;
    }

    if (chosenWindow == "hamming") {
        for (let n = 0; n < FRAMESIZE; n++) {
            chunk[n] = chunk[n] * (0.54 - 0.46 * Math.cos((2 * Math.PI * n) / (FRAMESIZE - 1)));
        }


    }
    if (chosenWindow == "blackman Harris") {
        for (let n = 0; n < FRAMESIZE; n++) {
            chunk[n] = chunk[n] * (0.35875 - 0.48829 * Math.cos((2 * Math.PI * n) / (FRAMESIZE - 1)) +
                0.14128 * Math.cos((4 * Math.PI * n) / (FRAMESIZE - 1)) -
                0.01168 * Math.cos((6 * Math.PI * n) / (FRAMESIZE - 1)));
        }
    }
    return chunk;

}

function addZeroes(frame) {
    const N = frame.length;
    if (N == nFFT) return frame;

    const numZeroes = nFFT - N;
    const leftZeroes = Math.floor(numZeroes / 2);

    const paddedFrame = new Float32Array(nFFT);

    paddedFrame.set(frame, leftZeroes);

    return paddedFrame;
}


function complexExp(angle) {
    return {
        real: Math.cos(angle),
        imag: Math.sin(angle)
    };
}



function createSpectrum(X) {

    ctx2D.clearRect(0, 0, canvas2D.width, canvas2D.height);
    //const windowRatio = nFFT / FRAMESIZE
    const K = X.length / 2; // Number of frequency bins
    const barWidth = canvas2D.width / K; // Width of each bar in the spectrum
    let freq = 0;
    // Normalize the magnitude values to fit the canvas height
    for (let i = 0; i < K / 2; i++) {
        //const freq = samplingFreq * (1 / K) * i; // Frequency (scaled for display purposes)

        const value = 9 * X[i] // Normalized magnitude


        // Use some color mapping for visualization
        const red = Math.min(255, value * 10);
        const green = Math.min(255, 255 - value * 5);
        const blue = Math.min(255, value * 2);
        ctx2D.fillStyle = `rgb(${red}, ${green}, ${blue})`;

        // Draw the bar
        ctx2D.fillRect(freq, canvas2D.height - value, barWidth, value);
        //ctx2D.fillRect(canvas2.width - barWidth, 0, barWidth, canvas2.height)
        freq += barWidth;
    }

}

function timeGraph(chunk) {


    const barWidth = 1; // Width of the bar
    const maxHeight = timeCanvas.height; // Centerline of the canvas
    // Create a copy of the current canvas
    const canvasCopy = ctxTime.getImageData(0, 0, timeCanvas.width, timeCanvas.height);

    // Clear the entire canvas
    ctxTime.clearRect(0, 0, timeCanvas.width, timeCanvas.height);

    // Redraw the copy, shifted left by 1 pixel
    ctxTime.putImageData(canvasCopy, -3, 0);

    // Clear the rightmost column
    ctxTime.clearRect(timeCanvas.width - (3 * barWidth), 0, 3 * barWidth, timeCanvas.height);
    for (let i = 3; i > 0; i--) {
        let value = maxHeight * chunk[i * 200]; // Supports positive and negative values

        // Draw the new bar based on the value's sign
        if (value >= 0) {
            // Positive bar: Draw above the centerline
            ctxTime.fillRect(
                timeCanvas.width - (i * barWidth), // X-position (rightmost edge)
                timeCanvas.height / 2 - value, // Y-position (centerline minus height)
                barWidth, // Width of the bar
                value // Positive height
            );
        } else {
            // Negative bar: Draw below the centerline
            ctxTime.fillRect(
                timeCanvas.width - (i * barWidth), // X-position (rightmost edge)
                timeCanvas.height / 2, // Y-position (centerline)
                barWidth, // Width of the bar
                -value // Negative height (make it positive for fillRect)
            );
        }
    }
}


function createMovingSpectrogram(X) {
    const barWidth = 1;
    const binHeight = canvasSpectrum.height / (nFFT / 2);
    ctxSpectrum.drawImage(canvasSpectrum, -barWidth, 0)
    ctxSpectrum.clearRect(canvasSpectrum.width - (barWidth), 0, barWidth, canvasSpectrum.height);
    //const windowRatio = nFFT / FRAMESIZE

    // Normalize the magnitude values to fit the canvas height
    const maxValue = Math.max(...X);
    const minValue = Math.min(...X);
    //const normalizedValues = X / maxMagnitude;
    //const normalizedValues = X.map(value => (value / maxMagnitude));
    //const maxNewMagnitude = Math.min(...normalizedValues);
    // Convert bin index to frequency
    const nyquist = SAMPLEFREQ / 2; // Nyquist frequency
    if (!melOn) {

        X.forEach((intensity, index) => {
            //(index / nFFT) * SAMPLEFREQ / 2
            const newIntensity = (intensity / SENS) - (CONTRAST);
            const frequency = (index / (nFFT / 2)) * nyquist;

            if (frequency <= 8000) {
                const yPosition = canvasSpectrum.height - (frequency / nyquist) * canvasSpectrum.height * SAMPLEFREQ / 16000;

                //const freqAxis = (index / K) * samplingFreq;
                //intensity = (intensity - 0.1) * 10
                ctxSpectrum.fillStyle = intensityToColor(newIntensity, maxValue, minValue);
                ctxSpectrum.fillRect(
                    (canvasSpectrum.width - barWidth),           // x-coordinate
                    yPosition,               // y-coordinate
                    barWidth,                        // width
                    binHeight                    // height
                );
            }
        });
    } else {
        const ratio = Math.ceil(SAMPLEFREQ / 16000);
        if (!melProcessed) {
            mel = melScale()
        }
        mel.forEach((intensity, index) => {
            const maxMel = 2595 * Math.log10(1 + (SAMPLEFREQ / 2) / 700); // Maximum mel value (Nyquist frequency)
            const val = (index / mel.length) * maxMel;
            const frequency = 700 * (10 ** ((val / 2595) - 1))
            if (intensity <= canvasSpectrum.height) {
                let melHeight = binHeight;
                if ((index > 0) && (index < nFFT / 2)) {
                    melHeight = (-mel[index - 1] + mel[index + 1]) / 2
                }

                ctxSpectrum.fillStyle = intensityToColor(X[index], maxValue, minValue);
                ctxSpectrum.fillRect(
                    canvasSpectrum.width - barWidth,           // x-coordinate
                    canvasSpectrum.height - mel[ratio * index],               // y-coordinate
                    barWidth,                        // width
                    melHeight * ratio                        // height
                );
            } else {
                //console.log(`Index: ${index}, Frequency: ${frequency.toFixed(2)} Hz`);
            }
        })





    }
}


function intensityToColor(intensity, maxValue, minValue) {
    const noiseThreshold = 0.1; // Define a threshold for noise (adjust as needed)
    const range = maxValue - minValue;
    let r, g, b;

    // Ensure the range is valid, avoid dividing by zero
    let normValue;
    if (chosenMagnitudeScale == "magnitude") {
        if (range == 0) {
            normValue = 0; // Handle edge case where max and min are equal
        } else if (intensity < noiseThreshold) {
            normValue = 0; // Treat very small values as noise
        } else {
            normValue = (intensity - minValue) / range;
            normValue = Math.max(0, Math.min(normValue, 1)); // Clamp to [0, 1]
        }
        const value = Math.floor((1 - normValue) * 255);


        if (chosenColourScheme == "greyScale") {
            // Map normalized intensity to grayscale
            r = value;
            g = value;
            b = value;
        } else if (chosenColourScheme == "neon") {
            if (range != 0) {
                r = 255 - value;
                g = value;
                b = 255;
            } else {
                r = 0;
                g = 255;
                b = 255;
            }
        } else if (chosenColourScheme == "heatedMetal") {
            // Map normalized intensity to a heated metal color scheme
            if (normValue <= 0.33) {
                // Red to Orange transition
                r = 255;
                g = Math.floor(255 * (normValue / 0.33));
                b = 0;
            } else if (normValue <= 0.66) {
                // Orange to Yellow transition
                r = 255;
                g = 255;
                b = Math.floor(255 * ((normValue - 0.33) / 0.33));
            } else {
                // Yellow to White transition
                r = 255;
                g = 255;
                b = 255;
            }
        }


    } else if (chosenMagnitudeScale == "deciBels") {
        const minIntensity = -150;
        const maxIntensity = 0;
        let normalized = Math.max(0, Math.min(1, (intensity - minIntensity) / (maxIntensity - minIntensity)));
        let normalizedPowered = Math.pow((normalized), 5)
        let value = Math.round((1 - normalizedPowered) * 255);
        // Map normalized value to grayscale
        // console.log(normalized)
        if (chosenColourScheme == "greyScale") {
            if (range != 0) {
                // Normalize intensity to a 0-255 scale (clamped to -150 dB to 0 dB)


                r = g = b = value;
            } else {
                r = g = b = 255; // Default black if range is 0
            }
        } else if (chosenColourScheme == "neon") {
            if (range != 0) {
                r = 255 - value;
                g = value;
                b = 255;
            } else {
                r = 0;
                g = 255;
                b = 255;
            }
        } else if (chosenColourScheme == "heatedMetal") {
            if (range != 0) {
                // Spread value across red to yellow to orange to red range
                r = Math.min(255, value * 2); // Red increases with value
                g = Math.min(255, value * 1.5); // Green is less intense than red
                b = Math.max(0, 255 - value); // Blue decreases as value increases
            } else {
                r = 255; // Default to full red when range is 0
                g = 255; // Yellow at full intensity
                b = 0;


            }
        }
    }


    return `rgb(${r}, ${g}, ${b})`;
}

function melScale() {
    const numBins = nFFT / 2;
    let maxFreq = 16000 / 2;
    //melOn = true;
    const maxMel = 2595 * Math.log10(1 + maxFreq / 700);
    let melBins = new Array(numBins);

    for (let i = 0; i < numBins; i++) {
        const freq = (i / numBins) * maxFreq;

        const currentMel = 2595 * Math.log10(1 + freq / 700);

        melBins[i] = (currentMel / maxMel) * canvasSpectrum.height;
    }
    melProcessed = true;


    return melBins;
}
/*
function intensityToColor(intensity, maxValue, minValue) {
    const normValue = (intensity - minValue) / (maxValue - minValue)
    const scaledValue = Math.floor((1 - normValue) * 255)

    const value = Math.floor((1 - intensity) * 255); // Map 0-1 to 0-255
    let r, g, b;


    if (chosenColourScheme == "greyScale") {

        r = scaledValue;
        g = scaledValue;
        b = scaledValue;

    }
    if (chosenColourScheme == "heatedMetal") {
        if (intensity <= 0.33) {
            // Red to Orange transition
            r = 255;
            g = Math.floor(255 * (intensity / 0.33)); // Scale green from 0 to 255
            b = 0;
        } else if (intensity <= 0.66) {
            // Orange to Yellow transition
            r = 255;
            g = 255;
            b = Math.floor(255 * (((1 - intensity) - 0.33) / 0.33)); // Scale blue from 0 to 255
        } else {
            // Yellow to White transition
            r = 255;
            g = 255;
            b = 255;
        }
    }

    return `rgb(${r}, ${g}, ${b})`;
}
*/

function drawAxisLabel() {
    const labelChunkeness = 2;
    const numLabels = 16;
    const labelWidth = 10;
    if (!melOn) {
        ctxAxis.clearRect(canvasAxis.width - 51, 0, 51, canvasAxis.height)

        const labelHeight = canvasAxis.height / numLabels;

        for (let n = 0; n <= numLabels; n++) {
            const label = `${((16000 / 2) / numLabels) * n} Hz`; // Example frequency labels
            ctxAxis.fillText(label, canvasAxis.width - 51, canvasAxis.height - n * labelHeight + 4); // Adjust position as needed

            ctxAxis.fillRect(
                canvasAxis.width - labelWidth,
                canvasAxis.height - n * labelHeight,
                labelWidth,
                labelChunkeness
            )

        }
    } else {
        ctxAxis.clearRect(canvasAxis.width - 51, 0, 51, canvasAxis.height)
        const melNum = Math.floor(mel.length / numLabels);

        for (let n = 0; n <= numLabels; n++) {
            const label = `${((16000 / 2) / numLabels) * n} Hz`; // Example frequency labels

            const melValue = mel[n * melNum];

            ctxAxis.fillText(label, canvasAxis.width - 51, canvasAxis.height - melValue); // Adjust position as needed
            ctxAxis.fillRect(
                canvasAxis.width - labelWidth,
                canvasAxis.height - melValue,
                labelWidth,
                labelChunkeness
            )
        }
    }
}