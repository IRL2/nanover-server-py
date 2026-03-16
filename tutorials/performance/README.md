FILES:
`\\Jeff\jeff-all\PROJECT FILES\PERFORMANCE TESTING`

Servers in `perf-demos.ipynb`
 * Choose sequence of input files (solvent cubes of various particle counts)
 * Run server
 * Run one of the recording cells (invisible or visible solvent)

Graphs in `perf-analysis.ipynb`

## Networking test
Requirements:
 * Install included `nanover-imd-vr-android.apk` on headset
 * Two WiFi connected machines, one for server, one for dummy clients

Running:
 * Run server from `perf-demos.ipynb` on first machine
 * Connect headset (use included `nanover-imd-vr-android.apk`)
 * Run n dummy clients (`dummy.py`) on second machine (you will see them as jumping headsets)
 * Choose output filename and run recording cell
 * Wear VR headset until server kicks you out

## Rendering test
Requirements:
 * Install included `nanover-imd-vr-android.apk` on headset
 * One WiFi connect machine for server

Running:
 * Run server from `perf-demos.ipynb`
 * Connect headset (use included `nanover-imd-vr-android.apk`)
 * Position box as desired
 * Choose output filename and run recording cell
 * Wear VR headset and maintain gaze until server kicks you out
