
   
 ## Aim  
 To perform Delta Modulation (DM) on a continuous signal and analyze its reconstruction.  
 ## 29/03/25
 ## Santhosh J
 ## 212223060248
 
 ## Tools Required  
 - Python (colab)
 - NumPy  
 - Matplotlib  
 - SciPy  
 
 ## Program  
 
 ### Delta Modulation and Demodulation  
 ```python
 import numpy as np
 import matplotlib.pyplot as plt
 from scipy.signal import butter, filtfilt
 
 # Parameters
 fs = 10000  # Sampling frequency
 f = 10  # Signal frequency
 T = 1  # Duration in seconds
 delta = 0.1  # Step size
 
 # Generate time vector
 t = np.arange(0, T, 1/fs)
 
 # Generate message signal (analog signal)
 message_signal = np.sin(2 * np.pi * f * t)  # Sine wave as input signal
 
 # Delta Modulation Encoding
 encoded_signal = []
 dm_output = [0]  # Initial value of the modulated signal
 prev_sample = 0
 
 for sample in message_signal:
     if sample > prev_sample:
         encoded_signal.append(1)
         dm_output.append(prev_sample + delta)
     else:
         encoded_signal.append(0)
         dm_output.append(prev_sample - delta)
     prev_sample = dm_output[-1]
 
 # Delta Demodulation (Reconstruction)
 demodulated_signal = [0]
 for bit in encoded_signal:
     if bit == 1:
         demodulated_signal.append(demodulated_signal[-1] + delta)
     else:
         demodulated_signal.append(demodulated_signal[-1] - delta)
 
 # Convert to numpy array
 demodulated_signal = np.array(demodulated_signal)
 
 # Apply a low-pass Butterworth filter
 def low_pass_filter(signal, cutoff_freq, fs, order=4):
     nyquist = 0.5 * fs
     normal_cutoff = cutoff_freq / nyquist
     b, a = butter(order, normal_cutoff, btype='low', analog=False)
     return filtfilt(b, a, signal)
 
 filtered_signal = low_pass_filter(demodulated_signal, cutoff_freq=20, fs=fs)
 
 # Plotting the Results
 plt.figure(figsize=(12, 6))
 
 plt.subplot(3, 1, 1)
 plt.plot(t, message_signal, label='Original Signal', linewidth=1)
 plt.title("Original Signal")
 plt.xlabel("Time [s]")
 plt.ylabel("Amplitude")
 plt.legend()
 plt.grid()
 
 plt.subplot(3, 1, 2)
 plt.step(t, dm_output[:-1], label='Delta Modulated Signal', where='mid')
 plt.title("Delta Modulated Signal")
 plt.xlabel("Time [s]")
 plt.ylabel("Amplitude")
 plt.legend()
 plt.grid()
 
 plt.subplot(3, 1, 3)
 plt.plot(t, filtered_signal[:-1], label='Demodulated & Filtered Signal', linestyle='dotted', linewidth=1, color='r')
 plt.title("Demodulated & Filtered Signal")
 plt.xlabel("Time [s]")
 plt.ylabel("Amplitude")
 plt.legend()
 plt.grid()
 
 plt.tight_layout()
 plt.show()
 ```
 
 ## Output Waveforms  
 - **Message Signal:** The original analog sine wave before Delta Modulation.
   
   
 ![Screenshot 2025-03-30 223816](https://github.com/user-attachments/assets/abcc4c42-e60c-48e3-ba88-2584b18e78aa)

 - **Delta Modulated Signal:** The stepwise encoded representation of the signal.
   
  
   ![Screenshot 2025-03-30 223833](https://github.com/user-attachments/assets/a34b17bf-8e5b-4880-82f4-ec42f7a8b521)

 
 - **Demodulated & Filtered Signal:** The reconstructed signal after filtering.
   
   
 ![Screenshot 2025-03-30 223859](https://github.com/user-attachments/assets/4d25cf74-3178-480a-848b-8f38cf612881)

 ## Results  
 - The analog message signal was successfully sampled and encoded using Delta Modulation. 
