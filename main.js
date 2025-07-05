const audioContext = new (window.AudioContext || window.webkitAudioContext)();


const SAMPLE_RATE = 16000;
const FRAME_SIZE = 1024;
const DEVICESAMPLERATE = audioContext.sampleRate;

const micWorker = new Worker('micWorker.js');
const fftWorker = new Worker('fftWorker.js');

worker.postMessage({ type: "config", sampleRate: SAMPLE_RATE,  deviceSampleRate: DEVICESAMPLERATE });
worker.postMessage({ type: "config", deviceSampleRate: DEVICESAMPLERATE });

let latestFFTData = [];

// Handle FFT output for visualization
fftWorker.onmessage = (e) => {
  latestFFTData = e.data;
};

// Create spectrogram visualization loop
function drawLoop() {
  requestAnimationFrame(drawLoop);
  if (latestFFTData.length > 0) {
    createMovingSpectrogram(latestFFTData, FRAME_SIZE);
  }
}

// Start drawing
drawLoop();

async function startAudioPipeline() {
  const audioContext = new AudioContext({ sampleRate: SAMPLE_RATE });

  // Load the AudioWorkletProcessor
  await audioContext.audioWorklet.addModule('micProcessor.js');

  // Create mic input
  const stream = await navigator.mediaDevices.getUserMedia({ audio: true });
  const source = audioContext.createMediaStreamSource(stream);

  // Create and connect audio worklet node
  const micNode = new AudioWorkletNode(audioContext, 'mic-processor');
  source.connect(micNode).connect(audioContext.destination); // destination is optional

  // Handle mic data
  micNode.port.onmessage = (e) => {
    const chunk = e.data;
    micWorker.postMessage(chunk, [chunk.buffer]);
  };

  // When mic worker finishes prepping audio frame
  micWorker.onmessage = (e) => {
    fftWorker.postMessage(e.data, [e.data.buffer]);
  };
}

startAudioPipeline().catch(console.error);
