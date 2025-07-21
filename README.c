/* READ ME -NOT UP TO DATE 20/07/2025
Created: 15/04/2025
Updated: 15/04/2025

Hugo Doughty
hdou141@aucklanduni.ac.nz

nFFT - the number of total fft points going into the fft, used in zero padding
audioBuffer - contains all the information for an audio wave

All the code for the real time spectrogram is stored in one file, index.js (not including 
the visual setup and appearance of the webpage, that is in index.html and index.css respectily). 
It runs off of two stages, running a pre recorded audio file, or accessing the microphone, 
displaying either in real time. 

---------------- AUDIO FILE INPUT -----------------------------------------
The run a file stage, essentially runs the audio file through a resampler function to account for 
the desired sampling frequency and sets up all the required audio buffer settings. Then it breaks 
the entire Audio file into chunks for the desired framesize and accounts for overlap and window type. 
To mimick real-time behaviour, it calculates the time for each chunk, and processes each one when required.
To process it, it zero pads the chunk and sends it through the foundation code of this software being the 
FFT function. This takes a chunk of samples, and returns the imaginary, real components of the 
frequency response and the raw magnitude and dB magnitude. These values are then sent to the 
createMovingSpectrogam function which handles the display of the spectrogram graph. On the side of this math 
part of the code, the audiofile is also played out loud, and if the real-time behaviour of the chunks is correct.
To the human eye, the spectrograph should display at approximatly the same time as the audio. 

---------------- MIRCROPHONE INPUT -----------------------------------------

The accessing microphone stage has a few more steps to set it up, but essentially mimicks the run a file stage.
The reasoning for extra steps, is becuase I wanted to code everything from scratch, maybe not the most efficient 
method, but I wanted to ensure complete control over all potential vaeriables, and this was the only way. 
To access the microphone, or an external device's audio, it needed to be stored in an array like a float array of 
32bits, this is raw values and not some classified audiobuffer with predefined audio charactersitcs, I had to 
assign these charactersitcs manually. Which also made it difficult as to create an audiobuffer and assign it to 
my analysers, meant I had to compete with the audiobuffer used for running a file. This is because I wanted to 
hijack or reuse the same processing code used to process an audiofile as, at the lowest level, the exact same 
process, just with a different input. 

The microphone stage grabs a number samples to create one chunk from the mircrophone, with the same length as 
desiredFrameSize, resamples and adds overlap. It then scales nFFT to the closest double of desiredFrameSize 
(This is to account for Nyquist) and then applies all the neccassary steps to preparing a chunk for FFT (apply 
window, add zeroes). The chunk is processed and the output is sent to createMovingSpectrogram. This process 
repeats infintly. 


Code Structure:
The stages are activated through interupts which trigger when a change occurs at the linked button 

----------------- AUDIO FILE INPUT ---------------------------------------
audioFileInput          - Is the File Input Interupt
setUpAudioBuffer()      - Sets audiobuffer correctly, audioBuffer contains audiofile info
executeFFTwithSync()    - Real-time fft application
    sliceIntoChunks()   
        addOverlap()
        applyWindow()
    processAudioFileChunks()     - Loops through and processes each chunk
        addZeroes()
        FFT()
        createMovingSpectogram()
        processAudioFileChunks()  --------- RECURSION

----------------- MICROPHONE INPUT ---------------------------------------

micButtonInput          - Is the interrupt
getMicData()            - Obtains MicData and sets everything up ready for proccessing
processMicInput()       - Recursice Loop to process and display mic data
    getMicData()
    ResampleBuffer()
    addOverlap()
    applyWindow()
    addZeroes()
    FFT()
    createMovingSpectrogram() 
    processMicInput()   -------------------- RECURSION 


The key difference between the two stages is that, Audio File Input slices the entire pre-existing audio file
into chunks and applies overlap, zeroes, and the window. Then it will process each one on the go and display them
real-time

Microphone access however, grabs each chunk as it comes and applies overlap (with previous chunk), zeroes and the 
window. Then will process the current chunk on the go and display it. 

For the microphone EVERYTHING is real-time. For audioFiles only the processing and display is real-time.




*/