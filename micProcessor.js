// micProcessor.js
class MicProcessor extends AudioWorkletProcessor {
  constructor() {
    super();
  }

  process(inputs, outputs, parameters) {
    const input = inputs[0]; // [channel][sample array]
    if (input.length > 0 && input[0].length > 0) {
    const buffer = new Float32Array(input[0]);
    this.port.postMessage(buffer, [buffer.buffer]);
    }
    return true; // Keep processor alive
  }
}

registerProcessor('mic-processor', MicProcessor);
