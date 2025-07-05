// micProcessor.js
class MicProcessor extends AudioWorkletProcessor {
  constructor() {
    super();
  }

  process(inputs, outputs, parameters) {
    const input = inputs[0]; // [channel][sample array]
    if (input.length > 0 && input[0].length > 0) {
      this.port.postMessage(input[0]); // Send raw mic samples to main thread
    }
    return true; // Keep processor alive
  }
}

registerProcessor('mic-processor', MicProcessor);
