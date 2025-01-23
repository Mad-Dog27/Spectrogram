const processAgainInput = document.getElementById('processAgain');
const micButtonInput = document.getElementById('Mic');

const melButtonInput = document.getElementById("melToggle")
const isMel = document.getElementById('isMel')

const recordButtonInput = document.getElementById('recordToggle')
const recordValue = document.getElementById('recordValue')

const playRecordButton = document.getElementById('playRecord')

const audioFileInput = document.getElementById('audioFileInput');
// Create an AudioContext
const audioContext = new (window.AudioContext || window.webkitAudioContext)();

/*
const FRAMESIZE = 2048; //time domain amount of samples taken
const nFFT = 2048; //frequency domain amount zeroes and values aquired through fft
const overlap = 1024;
*/
const FRAMESIZE = 2048 / 1; //time domain amount of samples taken
const nFFT = 4096 / 1; //frequency domain amount zeroes and values aquired through fft
const overlap = FRAMESIZE / 2;
const SPEED = 1;
let SCALE = 3;
let SENS = 1;
let CONTRAST = 0;
let recordOn = false;
let storedBuffer = [];
const container = document.getElementById("container");
const canvas = document.getElementById("canvas");
const canvas2 = document.getElementById("canvas2");

const canvasSpectrum = document.getElementById("canvasSpectrum");
const ctxSpectrum = canvasSpectrum.getContext("2d");


const scaleSlider = document.getElementById('scaleSlider');
const scaleSliderValue = document.getElementById('scaleSliderValue');

const sensSlider = document.getElementById('sensSlider');
const sensSliderValue = document.getElementById('sensSliderValue');

const contrastSlider = document.getElementById('contrastSlider');
const contrastSliderValue = document.getElementById('contrastSliderValue');


const windowSelect = document.getElementById('windowSelect');

const colourSchemeSelect = document.getElementById('colourSelect');

console.log(`Canvas dimensions: ${canvas2.width}x${canvas2.height}`);

console.log(`Canvas dimensions: ${canvasSpectrum.width}x${canvasSpectrum.height}`);




const ctx = canvas.getContext("2d");
const ctx2 = canvas2.getContext("2d");

// rectangular, hamming, blackman Harris
chosenWindow = "blackman Harris"
chosenColourScheme = 'greyScale'


let audioBuffer;
let fileUpload = null;
let filePlaying = null;
let finshedRT = null;
let micOn = null;

let mel = null;
let melOn = false;

//Adds an event listnener for the audioFileInput button, when the input file is changed the function will run after the file is changed.
audioFileInput.addEventListener('change', async (event) => {
    const file = event.target.files[0]; // file list of all the inputed files, [0] means only one file will be used
    if (!file) return;
    if (micOn) { console.log("Microphone is still recording"); return; }
    fileUpload = true;
    const arrayBuffer = await file.arrayBuffer(); //Read the file as an ArrayBuffer which is a binary representation of the audio file to use in the next line
    audioBuffer = await audioContext.decodeAudioData(arrayBuffer);//Uses the binary version to create an audio buffer

    melScale()
    // Process the audio buffer (e.g., generate a spectrogram)
    filePlaying = true;
    processAudioBuffer(audioBuffer);
});


processAgainInput.addEventListener('click', () => {
    //Event listener to process the audio file again, will happen on click
    if (!fileUpload) { console.log("No input file selected"); return; }
    if (filePlaying) { console.log("File is currently playing"); return; }
    if (micOn) { console.log("Mircophone Audio is currently playing"); return; }
    filePlaying = true;

    console.log("proccessing again")
    processAudioBuffer(audioBuffer);

})
melButtonInput.addEventListener('click', () => {
    //Event listener to process the audio file again, will happen on click

    if (melOn) {
        melOn = null;
        isMel.textContent = "off";
    } else {
        melOn = true;
        isMel.textContent = "on";
    }
    console.log("Mel toggled")

})

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

scaleSlider.addEventListener('input', () => {
    //Scale slider display and use, storign value in SCALE
    SCALE = 3 * scaleSlider.value / 50; //Converting to a smaller number for later maths
    scaleSliderValue.textContent = SCALE; // Update the display
    console.log(`Scale: ${SCALE}`);
});

sensSlider.addEventListener('input', () => {
    //Sensitivity slider display and use, on new input
    SENS = sensSlider.value; //Storing new value in SENS
    sensSliderValue.textContent = SENS; // Update the display
    console.log(`Senitivity: ${SENS}`);
});

contrastSlider.addEventListener('input', () => {
    //Contrast slider display and use, on new input
    CONTRAST = contrastSlider.value; //Storing each new value in CONTRAST
    contrastSliderValue.textContent = CONTRAST; //Updating display
    console.log(`Contrast: ${CONTRAST}`);
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

function processRecording() {
    audioBuffer = null;


    let combinedBuffer = (combineBuffers(storedBuffer))


    console.log(combinedBuffer)

    audioBuffer = audioContext.createBuffer(
        1, // Number of channels
        combinedBuffer.length, // Length of the buffer
        audioContext.sampleRate // Sample rate
    );

    // Copy the combined buffer into the AudioBuffer
    audioBuffer.getChannelData(0).set(combinedBuffer);

    // Create a source and play the buffer
    const source = audioContext.createBufferSource();
    source.buffer = audioBuffer;
    source.connect(audioContext.destination);
    source.start(0);

    processAudioBuffer(audioBuffer);    /*
    for (let currentChunk = 0; currentChunk < numChunk; currentChunk++) {
        const chunk = addZeroes(applyWindow(storedBuffer[currentChunk]))
        const result = fft(chunk);

        const dataMagnitude = result.map(bin => bin.magnitude);
        const slicedMagnitude = dataMagnitude.slice(0, nFFT / 2);
        createMovingSpectrogram(slicedMagnitude); //Create real-time spectrograph (MY CODE)

    }*/

}
function combineBuffers(buffers) {
    const totalLength = buffers.reduce((sum, buf) => sum + buf.length - overlap, 0) + overlap;
    const combinedBuffer = new Float32Array(totalLength);

    let offset = 0;
    buffers.forEach((buffer, index) => {
        for (let i = 0; i < buffer.length; i++) {
            combinedBuffer[offset + i] += buffer[i];
        }
        offset += buffer.length - overlap;
    });

    return combinedBuffer;
}



async function getMicData() {
    try {
        if (mel == null) { mel = melScale(); }
        // Requesting microphone acces and waiting, then creating a mediastream from it
        const stream = await navigator.mediaDevices.getUserMedia({ audio: true });
        const mediaStreamSource = audioContext.createMediaStreamSource(stream);

        // Create an AnalyserNode for real-time frequency domain analysis
        const analyser = audioContext.createAnalyser();
        analyser.fftSize = FRAMESIZE;
        mediaStreamSource.connect(analyser);

        // Array to store time-domain or frequency-domain data
        const dataArray = new Float32Array(analyser.frequencyBinCount);
        let chunkIndex = 0;
        // Function to process and visualize the microphone input 
        function processMicInput() {
            //New array/Buffer to store the audio samples
            let timeDomainBuffer = new Float32Array(FRAMESIZE);

            /*In the audio file processing, I store the audio data myself, but in this case I am not able to access
            the audio buffer, so am making use of getFloatTimeDomainData which stores the audio data from the analyser
            int my timeDomainBuffer. 
            NOTE: getFloatFrequencyData could also be used to obtain the frequency magnitudes But I prefer to use my maths*/
            analyser.getFloatTimeDomainData(timeDomainBuffer);

            const chunk = addZeroes(applyWindow(timeDomainBuffer));//Applying a window AND zero padding, the function above defaults to rectangular window
            const result = fft(chunk); //FAST FOURIER TRANSFORM

            const dataMagnitude = result.map(bin => bin.magnitude);
            const slicedMagnitude = dataMagnitude.slice(0, nFFT / 2)


            drawVisual(analyser); //Create first audio signal graph (NOT MY CODE)
            createSpectrum(dataMagnitude) //Create somewhat real-time spectrum (MY CODE)
            createMovingSpectrogram(slicedMagnitude); //Create real-time spectrograph (MY CODE)

            if (recordOn) {
                storedBuffer[chunkIndex] = timeDomainBuffer;
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


function processAudioBuffer(audioBuffer) {
    //Outputting the Audiobuffer characterestics
    console.log('Audio Buffer:', audioBuffer);
    console.log('Sample Rate:', audioBuffer.sampleRate);
    console.log('Number of Channels:', audioBuffer.numberOfChannels);
    console.log('Duration (s):', audioBuffer.duration);
    console.log(audioBuffer.length);
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

    let magnitudes = [];
    let result = [];
    let currentChunkIndex = 0;
    mel = melScale();
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

            // Fast Fourier transform of current chunk, cutting results in half then converting to magnitude
            result[currentChunkIndex] = fft(chunk);
            //result[currentChunkIndex] = result[currentChunkIndex].slice(0, FRAMESIZE / 2);

            //Updating all Graphs 
            const dataMagnitude = result[currentChunkIndex].map(bin => bin.magnitude);
            const slicedMagnitude = dataMagnitude.slice(0, nFFT / 2)

            const datadB = result[currentChunkIndex].map(bin => bin.dB);
            const slicedDb = datadB.slice(0, nFFT / 2)

            drawVisual(analyser)
            createSpectrum(slicedMagnitude);
            createMovingSpectrogram(slicedMagnitude);

            currentChunkIndex++; //incrementing chunk index
        }

        // Continue updating the spectrogram if there are more chunks to process
        if (currentChunkIndex < numChunks) {
            requestAnimationFrame(updateSpectrogram);
        } else {
            console.log("iedneifni"); filePlaying = null;
        }
    }

    // Start updating the spectrogram
    requestAnimationFrame(updateSpectrogram);
}

function sliceIntoChunks(audioBuffer) {

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
        combined[k + N / 2].dB = 20 * Math.log10(combined[k + N / 2].magnitude);
        //magnitude[k] = Math.sqrt(combined[k].real ** 2 + combined[k].imag ** 2);
    }

    return combined;
}
function complexExp(angle) {
    return {
        real: Math.cos(angle),
        imag: Math.sin(angle)
    };
}

function drawVisual(analyser) {
    const bufferLength = analyser.frequencyBinCount;//frequency bin count is a node property, bufferlength is the total num of frequency bins
    const dataArray = new Uint8Array(bufferLength);
    const barWidth = canvas.width / bufferLength; //divide the canvas width by freq bin amount to find the width of freq bin
    analyser.getByteFrequencyData(dataArray); //fills data array with freqeuncy data from audio signal, values represent loudness of each freq bin

    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx2.clearRect(0, 0, canvas.width, canvas.height);


    let x = 0;
    let x1 = canvas2.width / 2;
    let x2 = canvas2.width / 2;
    for (let i = 0; i < bufferLength; i++) {
        barHeight = dataArray[i]; //each frequency bin size

        const red = Math.min(255, barHeight * 10);
        const green = Math.min(255, 255 - barHeight);
        const blue = Math.min(255, barHeight * 2);
        ctx.fillStyle = `rgb(${red}, ${green}, ${blue})`;
        ctx.fillRect(x, (canvas.height - barHeight), barWidth, barHeight);//displaying freq bin size as a bar 
        /*
        ctx2.fillStyle = `rgb(${red}, ${green}, ${blue})`;
        ctx2.fillRect(x1, canvas.height - 5 * barHeight, barWidth, 5 * barHeight);//displaying freq bin size as a bar 
        ctx2.fillRect(x2, canvas.height - 5 * barHeight, barWidth, 5 * barHeight); // this will continue moving from left to right
        */
        x += barWidth;
        x1 -= barWidth; //moving x position to the right by bar width
        x2 += barWidth; //moving x position to the right by bar width

    }



}


function createSpectrum(X) {

    ctx2.clearRect(0, 0, canvas2.width, canvas2.height);
    //const windowRatio = nFFT / FRAMESIZE

    const K = X.length; // Number of frequency bins
    const barWidth = canvas2.width / K; // Width of each bar in the spectrum
    let freq = 0;
    // Normalize the magnitude values to fit the canvas height
    const maxMagnitude = Math.max(...X);
    const normalizedValues = X.map(value => (value / maxMagnitude) * canvas2.height);
    for (let i = 0; i < K / 2; i++) {
        //const freq = samplingFreq * (1 / K) * i; // Frequency (scaled for display purposes)

        const value = 8 * X[i] // Normalized magnitude


        // Use some color mapping for visualization
        const red = Math.min(255, value * 10);
        const green = Math.min(255, 255 - value * 5);
        const blue = Math.min(255, value * 2);
        ctx2.fillStyle = `rgb(${red}, ${green}, ${blue})`;

        // Draw the bar
        ctx2.fillRect(freq, canvas2.height - value, barWidth, value // Height of the bar
        );
        freq += barWidth;
    }

}

function createMovingSpectrogram(X) {
    const barWidth = 2;
    const binHeight = canvasSpectrum.height / (nFFT / 2);
    ctxSpectrum.drawImage(canvasSpectrum, -barWidth, 0)
    ctxSpectrum.clearRect(canvasSpectrum.width - barWidth, 0, barWidth, canvasSpectrum.height);
    //const windowRatio = nFFT / FRAMESIZE


    // Normalize the magnitude values to fit the canvas height
    const maxValue = Math.max(...X);
    const minValue = Math.min(...X);
    //const normalizedValues = X / maxMagnitude;
    //const normalizedValues = X.map(value => (value / maxMagnitude));
    //const maxNewMagnitude = Math.min(...normalizedValues);

    if (!melOn) {

        X.forEach((intensity, index) => {

            const newIntensity = (intensity / SENS) - (CONTRAST)

            //const freqAxis = (index / K) * samplingFreq;
            //intensity = (intensity - 0.1) * 10
            ctxSpectrum.fillStyle = intensityToColor(newIntensity, maxValue, minValue);
            ctxSpectrum.fillRect(
                canvasSpectrum.width - barWidth,           // x-coordinate
                canvasSpectrum.height - (index * binHeight * 2),               // y-coordinate
                barWidth,                        // width
                binHeight                    // height
            );

        });
    } else {
        mel.forEach((intensity, index) => {
            let melHeight = binHeight;
            if (index < nFFT - 1) {
                melHeight = mel[index + 1] - intensity
            }
            melHeight = melHeight / 1
            ctxSpectrum.fillStyle = intensityToColor(X[index], maxValue, minValue);
            ctxSpectrum.fillRect(
                canvasSpectrum.width - barWidth,           // x-coordinate
                canvasSpectrum.height - (2 * intensity),               // y-coordinate
                barWidth,                        // width
                melHeight                        // height
            );
        })


    }
}

function melScale() {
    const numBins = nFFT;
    let maxFreq = 8000;
    if (audioBuffer) {
        maxFreq = audioBuffer.sampleRate / 2;
    }
    const maxMel = 2595 * Math.log10(1 + maxFreq / 700);
    let melBins = new Array(numBins);

    for (let i = 0; i < numBins; i++) {
        const freq = (i / numBins) * maxFreq;

        const currentMel = 2595 * Math.log10(1 + freq / 700);

        melBins[i] = (currentMel / maxMel) * canvasSpectrum.height;
    }
    return melBins;
}

function intensityToColor(intensity, maxValue, minValue) {
    const normValue = (intensity - minValue) / (maxValue - minValue)
    const scaledValue = Math.floor((1 - normValue) * 255)

    const value = Math.floor((1 - intensity) * 255); // Map 0-1 to 0-255
    let r, g, b;

    /*if (chosenColourScheme == "greyScale") {
        if (maxValue > 0.1) {
            r = scaledValue;
            g = scaledValue;
            b = scaledValue;
        } else {
            r = 255;
            g = 255;
            b = 255;
        }
    }*/
    if (chosenColourScheme == "greyScale") {
        if (maxValue > 0) {
            r = value;
            g = value;
            b = value;
        } else {
            r = value;
            g = value;
            b = value;
        }
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




