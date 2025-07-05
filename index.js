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
              -Guassian Smoothing

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

NOTE: may crash your browser :) - Blame timeGraph() function - This is fixed 
*/
// WRITTEN ON THE 27/01/2025^^^






const optionWidth = "250px"
//Button INputs for inputing Data
const audioFileInput = document.getElementById('audioFileInput');
const processAgainInput = document.getElementById('processAgain');
processAgainInput.style.width = "120px"

const micButtonInput = document.getElementById('Mic');
micButtonInput.style.width = "120px"
const recordValue = document.getElementById('recordValue');
const playRecordButton = document.getElementById('playRecord');
//Button Inputs for altering graphs
const melButtonInput = document.getElementById("melToggle");
const isMel = document.getElementById('isMel');
const smoothButtonInput = document.getElementById('smoothButton');
const isSmooth = document.getElementById('isSmooth');
const recordButtonInput = document.getElementById('recordToggle');
const windowSelect = document.getElementById('windowSelect');
const colourSchemeSelect = document.getElementById('colourSelect');
const magnitudeSelect = document.getElementById('magnitudeSelect');
magnitudeSelect.style.width = optionWidth;
colourSchemeSelect.style.width = optionWidth;
windowSelect.style.width = optionWidth;
recordButtonInput.style.width = optionWidth;
smoothButtonInput.style.width = optionWidth;
melButtonInput.style.width = optionWidth;
playRecordButton.style.width = optionWidth;



const frameSizeSlider = document.getElementById('frameSizeSlider');
frameSizeSlider.style.width = optionWidth;
const frameSizeSliderValue = document.getElementById('frameSizeSliderValue');

const refSlider = document.getElementById("referenceSlider")
refSlider.style.width = "105px";
const refSliderValue = document.getElementById('referenceSliderValue')

const powerSlider = document.getElementById("powerSlider")
powerSlider.style.width = "105px";
const powerSliderValue = document.getElementById('powerSliderValue')

const contrastSlider = document.getElementById("contrastSlider")
const contrastSliderValue = document.getElementById('contrastSliderValue')


const widthSlider = document.getElementById('widthSlider');
widthSlider.style.width = optionWidth;
const widthSliderValue = document.getElementById('widthSliderValue');

const sampleFreqSlider = document.getElementById('sampleFreqSlider');
sampleFreqSlider.style.width = optionWidth;
const sampleFreqSliderValue = document.getElementById('sampleFreqSliderValue');

const overlapPercSlider = document.getElementById('overlapPercSlider');
overlapPercSlider.style.width = optionWidth;
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
let timeDiffs = []

ctxSpectrum.imageSmoothingEnabled = true;
//Global Constants
let FRAMESIZE = 128; //time domain amount of samples taken
let nFFT = 4096; //frequency domain amount zeroes and values aquired through fft
let overlapPercent = 0.25;
let overlap = Math.round(FRAMESIZE * overlapPercent);
const SPEED = 1;
let SAMPLEFREQ = 16000;
//Global Variables
let SCALE = 3;
let SENS = 1;
let CONTRAST = 1;
let recordOn = false;
let storedBuffer = [];
let noOverlapBuffer = [];
let audioBuffer;
let fileUpload = null;
let filePlaying = null;
let finshedRT = null;
let micOn = null;
let mel = null;
let melOn = false;
let smoothOn = null;
let melProcessed = false;
let isDB = null;
let REFERENCE = 2 ^ 15;
let WIDTH = 0.7;    //0.7
let HEIGHT = 0.7;  //0.49  DOESNT EFFECT - No point
let RecordProcessing = false;

let expectedChunkTime = FRAMESIZE / SAMPLEFREQ;
let frameUpdated = false;

let Size = 5;
let Sigma = 10.0;
let kernal = createGaussianKernel(Size, Sigma)

let shiftAccumulator = 0; // Global or function-scoped variable

let REF = 9.01;
let POW = 5;

let contrastOn = false;

let originalAudioBuffer;
chosenWindow = "blackman Harris"// rectangular, hamming, blackman Harris
chosenColourScheme = 'greyScale'
chosenMagnitudeScale = "magnitude"

canvasSpectrum.width = window.innerWidth * (WIDTH * 3) - 2;  // 70% of screen width minus borders
canvasSpectrum.height = window.innerHeight * (HEIGHT * 1) - 2;
canvasAxis.width = (canvasSpectrum.width + 2) / 3 - 2;
canvasAxis.height = canvasSpectrum.height;

timeCanvas.width = window.innerWidth * (WIDTH * 3) - 2;  // 70% of screen width minus borders
timeCanvas.height = window.innerHeight * 0.45 - 2;  // 45% of screen height minus borders

let timeCanvasRelativeWidth = (timeCanvas.width + 2) - 2
drawAxisLabel();

let filechunk = [];
let micchunk = [];

let intervalId = null;
let n =1;
let latestFrame = null; // or a ring buffer if you want to limit size
let effectiveChunkSize = FRAMESIZE - overlap;

let unregisteredUpdates = [];
let totalunreg = [];
let totalvals = [];


let toggle = false;
// ---------------------------vvv- INTERRUPTS -vvv----------------------------------------

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
    /*if (!fileUpload) { console.log("No input file selected"); return; }
    if (filePlaying) { console.log("File is currently playing"); return; }
    if (micOn) { console.log("Mircophone Audio is currently playing"); return; }
    filePlaying = true;
    let sum = 0;
    for (let i = 0; i < timeDiffs.length; i++) {
        sum += timeDiffs[i]
    }
    const avg = sum / timeDiffs.length
    console.log(avg)

    console.log("proccessing again")
    processAudioBuffer(audioBuffer);*/
    micWorker.terminate(); // if you're using one
     audioContext.suspend(); // to pause
// or
 audioContext.close(); // to fully stop everything
    toggle = true;
    fftWorker.terminate();

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
//mel scale input
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
//Actually for Gausian smoothing 
smoothButtonInput.addEventListener('click', () => {


    if (smoothOn) {
        smoothOn = null;
        isSmooth.textContent = "OFF"
    } else {
        smoothOn = true;
        isSmooth.textContent = "On"
    }
})
//magnitude select input
magnitudeSelect.addEventListener('change', (event) => {
    //Changing the window value, between Rectangular, hamming and blackman-harris
    chosenMagnitudeScale = event.target.value;
    console.log(`Window function selected: ${chosenMagnitudeScale}`);
    if (chosenMagnitudeScale == "deciBels") { isDB = true; }
    else {
        isDB = false;
    }
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
//frame size select
frameSizeSlider.addEventListener('input', () => {
    FRAMESIZE = parseInt(frameSizeSlider.value, 10);
    frameSizeSliderValue.textContent = FRAMESIZE;
    expectedChunkTime = FRAMESIZE / SAMPLEFREQ;
    console.log(expectedChunkTime)
    frameUpdated = true;
    intervalId = null;
    //nFFT = FRAMESIZE * 2; //frequency domain amount zeroes and values aquired through fft
    overlap = Math.round(FRAMESIZE * overlapPercent);
})
//reference input for dB
refSlider.addEventListener('input', () => {
    REF = refSlider.value;
    refSliderValue.textContent = REF;
})
// power input for dB 
powerSlider.addEventListener('input', () => {
    POW = powerSlider.value;
    powerSliderValue.textContent = POW;
})
//contrast slider input
contrastSlider.addEventListener('input', () => {
    CONTRAST = contrastSlider.value;
    contrastSliderValue.textContent = CONTRAST;
    contrastOn = true;

})

//width slider, for width of spectrogram displayed
widthSlider.addEventListener('input', () => {
    n = widthSlider.value;
    widthSliderValue.textContent = n; // Update the display
    console.log(`Width: ${n}`);
        intervalId = null;

    /*
    WIDTH = widthSlider.value;
    canvasSpectrum.width = window.innerWidth * (WIDTH*3) - 2;  // 70% of screen width minus borders

    timeCanvas.width = window.innerWidth * (WIDTH) - 2;  // 70% of screen width minus borders
    timeCanvasRelativeWidth = (timeCanvas.width + 2) * 3 - 2
    widthSliderValue.textContent = WIDTH; // Update the display
    console.log(`Width: ${WIDTH}`);
    */
});

//sample frequcy slider input
sampleFreqSlider.addEventListener('input', () => {//Function to update the INputed Sampling freq, this will improved freqeuency resolution but only untill you reach the original inputed frequency
    SAMPLEFREQ = sampleFreqSlider.value; //
    sampleFreqSliderValue.textContent = SAMPLEFREQ; // Update the display
    //adjust mel accordingly
    mel = melScale(); 
    //redraw the axis
    drawAxisLabel();
    expectedChunkTime = FRAMESIZE / SAMPLEFREQ;
    frameUpdated = true;
    intervalId = null;
    console.log(`FS: ${SAMPLEFREQ}`);
        console.log(expectedChunkTime)

});
// overlap input
overlapPercSlider.addEventListener('input', () => {
    overlapPercent = overlapPercSlider.value;
    overlapPercSliderValue.textContent = overlapPercent; // Update the display
    //recalculate amount of samples in overlap
    overlap = Math.round(FRAMESIZE * overlapPercent);
    console.log(overlap)

});

//Record toggle input
recordButtonInput.addEventListener('click', () => {
    if (recordOn) {
        recordOn = false;
        recordValue.textContent = "no"
    } else {
        storedBuffer = []; //reset recording
        recordOn = true;
        recordValue.textContent = "yes"
    }
})
// play the recording input
playRecordButton.addEventListener('click', () => {
    processRecording();
})

// -------------------------- INTERRUPTS DONE -----------------------------
const audioContext = new (window.AudioContext || window.webkitAudioContext)();


const SAMPLE__RATE = 16000;
const FRAME__SIZE = 1024;
const DEVICESAMPLERATE = audioContext.sampleRate;

const micWorker = new Worker('micWorker.js');
const fftWorker = new Worker('fftWorker.js');

micWorker.postMessage({ type: "config", sampleRate: SAMPLE__RATE,  deviceSampleRate: DEVICESAMPLERATE });
micWorker.postMessage({ type: "config", deviceSampleRate: DEVICESAMPLERATE });

let latestFFTData = new Float32Array(nFFT);
let fftResult = []
// Handle FFT output for visualization
fftWorker.onmessage = (e) => {
  const data = e.data;
    fftResult = []
    for (let i = 0; i < data.length; i += 2) {
    fftResult.push({ real: data[i], imag: data[i + 1] });
  }

    const computedResult = computeMagnitudeAndDB(fftResult, 1, false);

    const resultMagnitude = computedResult.map(bin => bin.magnitude);
    chosenValues = (resultMagnitude.slice(0, nFFT / 2));

    latestFFTData = chosenValues
    fftResult = 0
};

// Create spectrogram visualization loop
function drawLoop() {
    if (!toggle){
  requestAnimationFrame(drawLoop);
  if (latestFFTData.length > 0) {
    createMovingSpectrogram(latestFFTData, FRAMESIZE);
  }}
}

// Start drawing
drawLoop();

    async function startAudioPipeline() {
  const audioContext = new AudioContext({ sampleRate: SAMPLE__RATE });
    
  // Load the AudioWorkletProcessor
  await audioContext.audioWorklet.addModule('micProcessor.js');
    console.log("ifje")

  // Create mic input
  const stream = await navigator.mediaDevices.getUserMedia({ audio: true });
  const source = audioContext.createMediaStreamSource(stream);

  // Create and connect audio worklet node
  const micNode = new AudioWorkletNode(audioContext, 'mic-processor');
  source.connect(micNode).connect(audioContext.destination); // destination is optional
  
  // Handle mic data
  micNode.port.onmessage = (e) => {
  if (toggle) {
    micNode.disconnect();
    micNode.port.onmessage = null; 
    return; // 
  }

  const chunk = e.data;
  micWorker.postMessage(chunk);
};

  // When mic worker finishes prepping audio frame
  micWorker.onmessage = (e) => {
    fftWorker.postMessage(e.data, [e.data.buffer]);
  };
}

startAudioPipeline().catch(console.error);



let animationID = null
function renderSpectrogram() {
    if (micOn){
    animationID = requestAnimationFrame(renderSpectrogram);
    if (unregisteredUpdates !== null) {
        for (i = 0; i < unregisteredUpdates.length; i++) {

            createMovingSpectrogram(unregisteredUpdates[i], effectiveChunkSize);
            totalunreg.push([...unregisteredUpdates[i]])
         
        }
           unregisteredUpdates = []; // Optional: avoid redrawing same frame
    }
    } else {cancelAnimationFrame(animationID);
        console.log(totalunreg)
        console.log(totalvals)
        totalunreg = []
        totalvals = []
    }
}



//Grabs and processes the realtime mic data
async function getMicData() {
    /*Catherines IDEA: Create big chunk of raw audio data, then break into smaller desired chunks WHILE getting next big chunk, 
    process smaller chunks on the go. Maybe find another method of manually grabbing the timeBuffers instead of relying on Javascripts
    Methods/function. 
 
    If I could grab it manually I believe I would be able to accurattly procress and display the data. However in saying so, I have checked
    and the amount of data within each timeDomainBuffer is accurate, maybe its the rate at which they are being processed. The actual data IS
    different in magnitude compared to how it is in the Input File, this could be an issue.  
    */
    let ratio =  SAMPLEFREQ/  audioContext.sampleRate ; 
    let dataMagnitude, datadB;
    let prevChunkTime = 0;
    renderSpectrogram();

    try {
        if (mel == null) { }//mel = melScale(); }
        // Requesting microphone acces and waiting, then creating a mediastream from it
        const stream = await navigator.mediaDevices.getUserMedia({ audio: true });
        const mediaStreamSource = audioContext.createMediaStreamSource(stream);
        const deviceSampleRate = audioContext.sampleRate;
        //const ratio = SAMPLEFREQ / deviceSampleRate;
        let lowerPower = 1;
        let higherPower = 1;
        // Create an AnalyserNode for real-time frequency domain analysis
        const analyser = audioContext.createAnalyser();
        analyser.fftSize = 4096;
        mediaStreamSource.connect(analyser);
        let prevTimeDomainBuffer = [];
     
        let chunkIndex = 0;
        effectiveChunkSize = FRAMESIZE - overlap; //Subracting overlap, as that portion of chunk is accounted for in the next
        let closestFrameSize = FRAMESIZE;
        let neededFrameSize = FRAMESIZE / ratio;
        let closestNeededFrameSize = neededFrameSize;
        let runs = 0;
        let runTime = 0;
        const startTime = audioContext.currentTime;
        let chunkTime = 0;
        // Function to process and visualize the microphone input 
        function processMicInput() {
            
            // check if current time is larger then the expected chunk time
            if (chunkTime >= expectedChunkTime) {
                            //console.log(chunkTime)

            // calculate the closestneeded framesize (to the power of 2) to the users desired frame size,
            // this is used for resampling as, if lowereing Fs, more samples are needed to maintain same sized chunks
            if (frameUpdated) {
            while (lowerPower * 2 < neededFrameSize) { lowerPower <<= 1; }
            higherPower = lowerPower * 2;
            lowerPower >>= 1;
            closestNeededFrameSize = higherPower;
            higherPower = 1;
            lowerPower = 1;
            // calculating closest frame size to the power of 2 to the users desired. This is redundant right now
            while (lowerPower < FRAMESIZE) { lowerPower <<= 1; }
            higherPower = lowerPower;
            lowerPower >>= 1;
            closestFrameSize = higherPower;
            analyser.fftSize = closestNeededFrameSize; 
            frameUpdated = false;
            }

            if (closestFrameSize > nFFT) { 
                nFFT = closestFrameSize * 2; //frequency domain amount zeroes and values aquired through fft
            }
            //matching fftSize to the closestNeededFrameSize, this is becuase the analyser is only being used to collect the
            //time domain samples, there fftsize only effects the amount of samples taken. 

            

            //New array/Buffer to store the audio samples
            let timeDomainBuffer = new Float32Array(closestNeededFrameSize); //Larger than expected, see insied resample
            
            /*In the audio file processing, I store the audio data myself, but in this case I am not able to access
            the audio buffer, so am making use of getFloatTimeDomainData which stores the audio data from the analyser
            int my timeDomainBuffer. 
            NOTE: getFloatFrequencyData could also be used to obtain the frequency magnitudes But I prefer to use my maths*/
            analyser.getFloatTimeDomainData(timeDomainBuffer);

            const resampledTimeDomainBuffer = resampleMicBuffer( //resampling
                timeDomainBuffer,
                deviceSampleRate,
                SAMPLEFREQ,
                closestFrameSize
            )
            let newTimeDomainBuffer = new Float32Array(FRAMESIZE) //adding overlap to the resampled buffer
            newTimeDomainBuffer = addOverLap(resampledTimeDomainBuffer, prevTimeDomainBuffer);
            prevTimeDomainBuffer = resampledTimeDomainBuffer;

            const chunk = addZeroes(applyWindow(newTimeDomainBuffer, FRAMESIZE));//Applying a window AND zero padding, the function above defaults to rectangular window

            const result = fft(chunk); //FAST FOURIER TRANSFORM
            
            runs++;

            //magnitude selection 
            if (chosenMagnitudeScale == "magnitude") {
                const data = computeMagnitudeAndDB(result, 1, false);

                dataMagnitude = data.map(bin => bin.magnitude);
                chosenValues = (dataMagnitude.slice(0, nFFT / 2));
                
            } else {
                const data = computeMagnitudeAndDB(result, 1, true);

                datadB = data.map(bin => bin.dB);
                console.log(datadB)

                chosenValues = datadB.slice(0, nFFT / 2)
            }
                                    let time1 = audioContext.currentTime;

            unregisteredUpdates.push([...chosenValues]); // push a copy to avoid mutation
            let time2 = audioContext.currentTime +1;
            console.log("time ", (time2 - time1))
            //createMovingSpectrogram(chosenValues)
            //createSpectrum(dataMagnitude);
            //createMovingSpectrogram(chosenValues, effectiveChunkSize);
            if (recordOn) { //record mic input 
                if (!storedBuffer[chunkIndex]) {
                    noOverlapBuffer[chunkIndex] = resampledTimeDomainBuffer;
                    storedBuffer[chunkIndex] = newTimeDomainBuffer;
                } else {
                    console.warn("Chunk index already filled, potential overwrite detected.");
                }
                chunkIndex++;
            }
            }

            runTime = audioContext.currentTime - startTime;
            chunkTime = runTime - prevChunkTime; //calculating how long this chunk has occured.
            prevChunkTime += chunkTime; //sum of all prev chunk times

            if (intervalId == null) {
             intervalId = setInterval(() => {
            if (micOn) {
                processMicInput();
            } else {
                
                clearInterval(intervalId);
                intervalId = null;
            mediaStreamSource.disconnect();
                analyser.disconnect();
                stream.getTracks().forEach(track => track.stop());
                console.log('Stream processing stopped.');
                console.log("runTime: ", runTime);
                console.log("runs: ", runs);
                console.log("average chunk time: ", runTime / runs);
                console.log("expected chunk time: ", expectedChunkTime);

            }
            }, expectedChunkTime * 1000 * n); // e.g. 32ms = 0.032 * 1000
            }

            
           
     


        }

        // Start processing
        processMicInput();
    } catch (error) {

        console.error("Error accessing microphone:", error); //Might be better then using console.log
    }
}

// addoverlap to a buffer using the previous buffer
function addOverLap(timeDomainBuffer, prevTimeDomainBuffer) { //WRONG OVERLAP SHOUDLNT INCREASE FRAMESIZE
    
    const prevLength = prevTimeDomainBuffer.length;
    //new overlap in terms of samples
    const newOverlap = Math.floor(overlapPercent * FRAMESIZE);
    
    let newCurrentBuffer = new Float32Array(FRAMESIZE)
    if (prevLength == 0) { // if no previous buffer, then no overlap applied 
        newCurrentBuffer.set(timeDomainBuffer.subarray(0, newCurrentBuffer.length - newOverlap), newOverlap);
        return newCurrentBuffer;
    }

    for (let i = 0; i < newOverlap; i++) { // add end overlap amount of samples of previous buffer to the start of the new buffer
        newCurrentBuffer[i] = prevTimeDomainBuffer[prevLength - newOverlap + i]
    }
    // add the rest of the buffer
    newCurrentBuffer.set(timeDomainBuffer.subarray(0, newCurrentBuffer.length - newOverlap), newOverlap)

    return newCurrentBuffer;
}

// resampling an audiobuffer for File Input only
async function resampleAudio(audioBuffer) {
    const originalSampleFreq = audioBuffer.sampleRate;
    // if no resampling needed
    if (originalSampleFreq == SAMPLEFREQ) {
        return audioBuffer;
    }

    // new audioconext to do the resampling, offline means its local and not conflicting with global audio context
    const offlineContext = new OfflineAudioContext(
        audioBuffer.numberOfChannels,
        Math.ceil(audioBuffer.length * (SAMPLEFREQ / originalSampleFreq)),
        SAMPLEFREQ
    );
    // source and audio buffer attached to that offline audiocontext
    const source = offlineContext.createBufferSource();
    source.buffer = audioBuffer;
    source.connect(offlineContext.destination);

    source.start(0);

    const newBuffer = await offlineContext.startRendering();
    return newBuffer; //new resampled buffer

}

//resampling for mic buffer, (it cant use audioconexts)
function resampleMicBuffer(buffer, originalRate, targetRate, closestFrameSize) {//This works
    //Buffer is going to be ratio times larger then need be, this is becuase of interpolation, lower freqs skip a certain
    // amount of samples, but still needs to be same size, therefore original bufffer needs to be larger

    const ratio = originalRate / targetRate; // deviceFS / sampleFREQ (user chosen)
    let newBuffer = new Float32Array(FRAMESIZE);
    if (originalRate == targetRate) { return buffer; } // if no resampling needed
    for (let i = 0; i < FRAMESIZE; i++) { //newLength should be buffers original size(desired size)
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

// Process recording from microphone - still in testing
function processRecording() {
    for (let chunkIndex = 0; chunkIndex < storedBuffer.length; chunkIndex++) {
        const chunk = addZeroes(applyWindow(storedBuffer[chunkIndex], FRAMESIZE));//Applying a window AND zero padding, the function above defaults to rectangular window
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
    }
}
   
// process file input audio buffer (setup)
async function processAudioBuffer(audioBuffer) {
    //Outputting the Audiobuffer characterestics
    console.log('Audio Buffer:', audioBuffer);
    console.log('Sample Rate:', audioBuffer.sampleRate);
    console.log('Number of Channels:', audioBuffer.numberOfChannels);
    console.log('Duration (s):', audioBuffer.duration);
    console.log(audioBuffer.length);
    originalAudioBuffer = audioBuffer;
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



    /*For the next step, the audioBuffer needs to be cut into respective chunks for audioProccessing, 
    first the chunk characteristics need to be determined*/
    const sampleRate = audioBuffer.sampleRate;
    console.log(sampleRate)
    const effectiveChunkSize = FRAMESIZE - overlap; //Subracting overlap, as that portion of chunk is accounted for in the next
    const chunkDuration = effectiveChunkSize / SAMPLEFREQ; // Duration of one chunk in seconds
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
    let prevTime = startTime;
    console.log(startTime)
    let n = 0
    let m = 0;
    let alter = true;
    let runs = 0;
    let runTime = 0;
    // Function to update the spectrogram
    function updateSpectrogram() {
        runs++
        // Calculate the expected chunk index based on current playback time
        const elapsedTime = audioContext.currentTime - startTime;
        const expectedChunkIndex = Math.floor(elapsedTime / chunkDuration);
        //EACH ITERATION IS APPROXX 0.018s
        if (n != 0) {
            timeDiffs[n] = (audioContext.currentTime - startTime - prevTime)
        } else {
            timeDiffs[0] = 0;
        }
        n++
        // Process and render all chunks up to the expected index -- ADD A TIME BASED CONSTRAINT
        while (currentChunkIndex <= expectedChunkIndex && currentChunkIndex < numChunks) { //If chunk processing is slower then expected and chunk index has not exceded total number of chunks
            const chunk = addZeroes(chunks[currentChunkIndex]); //Zero Padding
            m++

            // Fast Fourier transform of current chunk, cutting results in half then converting to magnitude
            result = fft(chunk);

            // magnitude select
            if (!isDB) {
                //const dataMagnitude = result.map(bin => bin.magnitude);
                //chosenValues = (dataMagnitude.slice(0, nFFT / 2));
                const data = computeMagnitudeAndDB(result);
                const dataMagnitude = data.map(bin => bin.magnitude);
                chosenValues = (dataMagnitude.slice(0, nFFT / 2));

            } else {
                const datadB = result.map(bin => bin.dB);
                chosenValues = datadB.slice(0, nFFT / 2)
            }
            /*
            if (b == 20) {
                console.log(chunk)
                filechunk = chunk;
            }
            b++;*/
            if (smoothOn) { // NOT TIME ON, gaussian smoothing
                const smoothChosenValues = applyGaussianFilter(chosenValues, kernal)

                //createSpectrum(chosenValues);
                createMovingSpectrogram(smoothChosenValues, effectiveChunkSize);

            } else { // normal display
                //createSpectrum(chosenValues);
                createMovingSpectrogram(chosenValues, effectiveChunkSize);
            }
            currentChunkIndex++; //incrementing chunk index
            
        }

        // Continue updating the spectrogram if there are more chunks to process
        if (currentChunkIndex < numChunks) {
            requestAnimationFrame(updateSpectrogram);
        } else {
            filePlaying = null;
            console.log("runTime: ", runTime);
            console.log("runs: ", runs);

        }

        runTime = audioContext.currentTime - startTime;

    }

    // Start updating the spectrogram

    requestAnimationFrame(updateSpectrogram);

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
function normalise(values) { // this function works perfectly, BUT commented out because normalisation doesnt work when noise is only presented
    return values;

    maxVal = Math.max(...values);
    //if (maxVal < 0.008) return values;
    let normValues = new Array(values.length)
    console.log(maxVal)
    for (let i = 0; i < values.length; i++){
        normValues[i] = values[i]/maxVal;
    }
    return normValues;
} 
function createGaussianKernel(size, sigma) { // gaussian smoothing CAN BE ADJUSTED
    const newKernel = new Float32Array(size);
    const mean = Math.floor(size / 2);
    let sum = 0; // For normalization

    for (let x = 0; x < size; x++) {
        const exponent = -((x - mean) ** 2) / (2 * sigma ** 2);
        newKernel[x] = Math.exp(exponent);
        sum += newKernel[x];
    }

    // Normalize the kernel so the sum equals 1
    for (let x = 0; x < size; x++) {
        newKernel[x] /= sum;
    }

    return newKernel;
}
function applyGaussianFilter(chosenValues, kernal) { // Applying gaussian smoothing
    const smoothedValues = new Float32Array(chosenValues.length);
    const half = Math.floor(kernal.length / 2);

    for (let i = 0; i < chosenValues.length; i++) {
        let sum = 0;

        for (let j = -half; j <= half; j++) {
            const index = i + j;

            // Ensure the index is within bounds
            if (index >= 0 && index < chosenValues.length) {
                sum += chosenValues[index] * kernal[half + j];
            }
        }

        smoothedValues[i] = sum;
    }

    return smoothedValues;
}
/*
function fft(input, start = 0, stride = 1, N = input.length) {
        if ((N & (N - 1)) !== 0) {
            console.error("Input length must be a power of 2");
            return [];
        }
    
        // Base case
        if (N === 1) {
            return [{ real: input[start], imag: 0 }];
        }
    
        // Recursive FFT on even and odd indices
        const even = fft(input, start, stride * 2, N / 2);
        const odd = fft(input, start + stride, stride * 2, N / 2);
    
        const combined = Array(N).fill(0).map(() => ({ real: 0, imag: 0 }));
    
        for (let k = 0; k < N / 2; k++) {
            const angle = (-2 * Math.PI * k) / N;
            const twiddle = {
                real: Math.cos(angle),
                imag: Math.sin(angle)
            };
    
            const t = {
                real: twiddle.real * odd[k].real - twiddle.imag * odd[k].imag,
                imag: twiddle.real * odd[k].imag + twiddle.imag * odd[k].real
            };
    
            combined[k] = {
                real: even[k].real + t.real,
                imag: even[k].imag + t.imag
            };
            combined[k + N / 2] = {
                real: even[k].real - t.real,
                imag: even[k].imag - t.imag
            };
        }
    
        return combined;
    }*/
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

function sliceIntoChunks(audioBuffer) {
    let tempOverlap = overlap;

    if (RecordProcessing) {//DOUBLE CHECK THIS
        overlap = 0;
    }
    //averaging the audio if there are two channels
    const monoChannel = audioBuffer.getChannelData(0);
    if (audioBuffer.numberOfChannels >= 2) {
        const samples1 = audioBuffer.getChannelData(1);
        for (let j = 0; j < audioBuffer.length; j++) {
            monoChannel[j] = (monoChannel[j] + samples1[j]) / 2
        }
    }
    const numChunks = Math.floor((audioBuffer.length - FRAMESIZE) / (FRAMESIZE - overlap)) + 1;
    const chunks = [];

    //storing each chunk into array of chunks, finding the start and end index each iteration to do so
    for (let i = 0; i < numChunks; i++) {
        const start = i * (FRAMESIZE - overlap);
        const realStart = parseInt(start, 10);

        const end = parseInt(realStart + FRAMESIZE, 10);
        const chunk = monoChannel.slice(realStart, end);

        const windowedChunk = applyWindow(chunk, FRAMESIZE);
        //const chunk = audioBuffer.getChannelData(0).slice(offset, offset + bufferLength);
        chunks.push(windowedChunk);
    }
    overlap = tempOverlap; // cant remeber why this is needed 
    return chunks
}
// apply a window function to each chunk BEFORE zero padding
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

function addZeroes(frame) { // add zeroes to match nFFT value 
    let N = frame.length;
    if (N > nFFT) N = nFFT; // edge case, SHOULD NEVER HAPPEN

    if (N == nFFT) return frame; // no change needed
    const numZeroes = nFFT - N;
    const leftZeroes = Math.floor(numZeroes / 2); // zero padding to left of samples

    const paddedFrame = new Float32Array(nFFT);

    paddedFrame.set(frame, leftZeroes); // right zeropadding automatically done when creating new array ^^

    return paddedFrame;
}

// calculating real and imag based off angle
function complexExp(angle) {
    return {
        real: Math.cos(angle),
        imag: Math.sin(angle)
    };
}


// realtime graph to show which freqs are loud - just for visual not needed
function createSpectrum(X) {

    ctx2D.clearRect(0, 0, canvas2D.width, canvas2D.height);
    //const windowRatio = nFFT / FRAMESIZE
    const K = X.length; // Number of frequency bins
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
        ctx2D.fillRect(2 * freq, canvas2D.height - value, 2 * barWidth, value);
        //ctx2D.fillRect(canvas2.width - barWidth, 0, barWidth, canvas2.height)
        freq += barWidth;
    }

}


// function to display the audio wave - VERY LAGGY
function timeGraph(chunk) {
    const num = 3;
    const frameRatio = Math.floor(FRAMESIZE / num)
    const barWidth = 1; // Width of the bar
    const maxHeight = timeCanvas.height; // Centerline of the canvas
    // Create a copy of the current canvas
    const canvasCopy = ctxTime.getImageData(0, 0, timeCanvas.width, timeCanvas.height);

    // Clear the entire canvas
    ctxTime.clearRect(0, 0, timeCanvas.width, timeCanvas.height);

    // Redraw the copy, shifted left by 1 pixel
    ctxTime.putImageData(canvasCopy, -num, 0);
    // Clear the rightmost column
    ctxTime.clearRect(timeCanvas.width - (num * barWidth), 0, num * barWidth, timeCanvas.height);
    for (let i = num; i > 0; i--) {

        let value = maxHeight * chunk[i * frameRatio]; // Supports positive and negative values
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

// SPECTROGRAM DISPLAY FUNCTION
/*function createMovingSpectrogram(X, effectiveChunkSize) {//Mic is having scaling issues because sample freq has no effect on it
    let ratio = SAMPLEFREQ / 16000; // baseline of fs 16000, need a reference, however I am not too sure I like this - CHECK
    let barWidth = 1;
 
    ctxSpectrum.drawImage(canvasSpectrum, -barWidth, 0); // shift graph

    const binHeight = canvasSpectrum.height / (nFFT / 2);
    ctxSpectrum.clearRect(canvasSpectrum.width - (barWidth), 0, barWidth, canvasSpectrum.height); // delete first column

    const maxValue = Math.max(...X);
    const minValue = Math.min(...X);

    const nyquist = SAMPLEFREQ / 2; // Nyquist frequency
    if (!melOn) {
        X.forEach((intensity, index) => {
            const newIntensity = intensity; // no change
            const frequency = (index / (nFFT / 2)) * nyquist; //caluclating height of specific rectangle
                const yPosition = canvasSpectrum.height - (frequency / nyquist) * canvasSpectrum.height * ratio;

                ctxSpectrum.fillStyle = intensityToColor(newIntensity, maxValue, minValue); // calculate the colour
                ctxSpectrum.fillRect( // add new column 
                    (canvasSpectrum.width - barWidth),           // x-coordinate
                    yPosition,               // y-coordinate
                    barWidth,                        // width
                    binHeight * ratio                   // height
                );
            
        });
    } else { // if mel scaling 
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
                    melHeight = (-mel[index - 1] + mel[index + 1])
                }

                ctxSpectrum.fillStyle = intensityToColor(X[index], maxValue, minValue);
                ctxSpectrum.fillRect(
                    canvasSpectrum.width - barWidth,           // x-coordinate
                    canvasSpectrum.height - mel[ratio * index],               // y-coordinate
                    barWidth,                        // width
                    melHeight                // height
                );
            } else {
                //console.log(`Index: ${index}, Frequency: ${frequency.toFixed(2)} Hz`);
            }
        })

}
}*/
function createMovingSpectrogram(X, effectiveChunkSize) {
    const fs = SAMPLEFREQ;
    const barWidth = 2;
    const height = canvasSpectrum.height;
    const width = canvasSpectrum.width;
    const nyquist = fs / 2;
    const binHeight = height / (nFFT / 2);
    const ratio = fs / 16000;

    // Shift canvas to the left
    ctxSpectrum.drawImage(canvasSpectrum, -barWidth, 0);
    ctxSpectrum.clearRect(width - barWidth, 0, barWidth, height);

    // Precompute max/min values once (could be passed in from FFT too)
    let maxValue = -Infinity;
    let minValue = Infinity;
    for (let i = 0; i < X.length; i++) {
        const val = X[i];
        if (val > maxValue) maxValue = val;
        if (val < minValue) minValue = val;
    }

    if (!melOn) {
        // Linear frequency scale
        const xCoord = width - barWidth;

        for (let i = 0; i < X.length; i++) {
            const intensity = X[i];
            const y = height - (i / (nFFT / 2)) * height * ratio;
            const h = binHeight * ratio;

            ctxSpectrum.fillStyle = intensityToColor(intensity, maxValue, minValue);
            ctxSpectrum.fillRect(xCoord, y, barWidth, h);
        }
    } else {
        // MEL frequency scale
        if (!melProcessed) {
            mel = melScale();
        }

        const maxMel = 2595 * Math.log10(1 + nyquist / 700);
        const xCoord = width - barWidth;

        for (let i = 1; i < mel.length - 1; i++) {
            const melY = height - mel[ratio * i];
            const melHeight = mel[i + 1] - mel[i - 1]; // Local bandwidth approximation

            ctxSpectrum.fillStyle = intensityToColor(X[i], maxValue, minValue);
            ctxSpectrum.fillRect(xCoord, melY, barWidth, melHeight);
        }
    }
}


// calculates the colour based off value magnitude -  Not a fan of how it calculates normVal
function intensityToColor(intensity, maxValue, minValue) {
    const noiseThreshold = 0.01; // Define a threshold for noise (adjust as needed)
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


        const value = Math.floor((1 - intensity) * 255); // Map 0-1 to 0-255
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
        } else if (chosenColourScheme == "fancy") {

            if (normValue <= 0.2) {
                // Dark red to red (low magnitude)
                r = 40 + Math.floor(128 * (normValue / 0.2)); // 0 to 128
                g = 0;
                b = 0;
            } else if (normValue <= 0.3) {
                r = 90 + Math.floor(128 * (normValue / 0.2)); // 0 to 128
                g = 0;
                b = 0;
            } else if (normValue <= 0.4) {
                // Red to yellow
                r = 158 + Math.floor(127 * ((normValue - 0.2) / 0.2)); // 128 to 255
                g = Math.floor(255 * ((normValue - 0.2) / 0.2)); // 0 to 255
                b = 0;
            } else if (normValue <= 0.5) {
                // Yellow to cyan
                r = Math.floor(255 * (1 - (normValue - 0.4) / 0.2)); // 255 to 0
                g = 255;
                b = Math.floor(255 * ((normValue - 0.4) / 0.2)); // 0 to 255
            } else if (normValue <= 0.8) {
                // Cyan to blue
                r = 0;
                g = Math.floor(255 * (1 - (normValue - 0.6) / 0.2)); // 255 to 0
                b = 255;
            } else {
                // Blue to dark blue (high magnitude)
                r = 0;
                g = 0;
                b = Math.floor(255 * (1 - (normValue - 0.8) / 0.2)); // 255 to 0
            }

        }


    } else if (chosenMagnitudeScale == "deciBels") { // if decibels 
        const minIntensity = -150;
        const maxIntensity = 0;
        let normalized = Math.max(0, Math.min(1, (intensity - minIntensity) / (maxIntensity - minIntensity)));
        let normalizedPowered = Math.pow((normalized), POW) //Questionable alteration of Decibell scale, could change min and max
        let value = Math.round((1 - normalizedPowered) * 255);

        if (chosenColourScheme == "greyScale") {
            if (contrastOn) {
                /*if (value < 100) {
                    val = 100 - value
                    value -= val * CONTRAST;

                } */
                if (value > 150) {
                    val = value - 150
                    value = value + (value - 150) * CONTRAST;
                }

            }

            if (range != 0) {

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

// mel scaling, only needs to occur once every sample freq change - calulates specfic index etc 
function melScale() {
    const numBins = nFFT / 2;
    let maxFreq = SAMPLEFREQ / 2;
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

// draw the freq axis labels, before anything happens, 
function drawAxisLabel() {
    const labelChunkeness = 2;
    const numLabels = 20;
    const labelWidth = 10;
    if (!melOn) { // if normal
        ctxAxis.clearRect(canvasAxis.width - 51, 0, 51, canvasAxis.height)
        const labelHeight = canvasAxis.height / numLabels;
        console.log(canvasAxis.width)
        for (let n = 1; n <= numLabels; n++) {
            const label = `${((16000 / 2) / numLabels) * n} Hz`; // Example frequency labels
            ctxAxis.fillText(label, canvasAxis.width - 51, canvasAxis.height - n * labelHeight + 4); // Adjust position as needed

            ctxAxis.fillRect(
                canvasAxis.width - labelWidth,
                canvasAxis.height - n * labelHeight,
                labelWidth,
                labelChunkeness
            )

        }
    } else { // if mel
        ctxAxis.clearRect(canvasAxis.width - 51, 0, 51, canvasAxis.height)
        const melNum = Math.floor(mel.length / numLabels);

        for (let n = 0; n <= numLabels; n++) {
            const label = `${((16000 / 2) / numLabels) * n} Hz`; // Example frequency labels

            const melValue = mel[n * melNum];

            ctxAxis.fillText(label, canvasAxis.width - 51, canvasAxis.height - melValue + 4); // Adjust position as needed
            ctxAxis.fillRect(
                canvasAxis.width - labelWidth,
                canvasAxis.height - melValue,
                labelWidth,
                labelChunkeness
            )
        }
    }
}