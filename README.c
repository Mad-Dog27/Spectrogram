/* READ ME

All the code for the real time spectrogram is stored in one file, index.js (not including 
the visual setup and appearance of the webpage). It runs off of two stages, running 
a pre recorded audio file, or accessing the microphone, displaying either in real time. 



The run a file stage, essentially run the audio file through a resampler function to account for 
the desired sampling frequency and sets up all the required audio buffer settings. Then it breaks 
the entire Audio file into chunks for the desired framesize and accounts for overlap and window type. 
To mimick real-time behaviour, it calculates the time for each chunk, and processes each one when required.
To process it, it zero pads the chunk and sends it through the foundation code of this software being the 
FFT function. This takes a chunk of samples, and returns the imaginary, real components of the 
frequency response and the raw magnitude and dB magnitude. These values are then sent to the 
createMovingSpectrogam function which handles the display of the spectrogram graph. On the side of this math 
part of the code, the audiofile is also played out loud, and if the real-time behaviour of the chunks is correct,
to the human eye, the spectrograph should display at the approximatly same time as the audio. 


The accessing microphone stage has a few more steps to set it up, but essentially mimicks the run a file stage.
The reasoning for extra steps, is becuase I wanted to code everything from scratch, maybe not the most efficient 
method, but I wanted to ensure complete control over all potential vaeriables, and this was the only way. 
To access the microphone, or an external device's audio, it needed to be stored in an array like a float array of 
32bits, this is raw values and not some classified audiobuffer with predefined audio charactersitcs, I had to 
assign these charactersitcs manually. Which also made it difficult as to create an audibuffer and assign it to 
my analysers, meant I had to compete with the audiobuffer used for running a file. This is because I wanted to 
hijack or reuse the same processing code used to process an audiofile as, at the lowest level, the exact same 
process, just with a different input. 


*/