## Settings files for trap frequency Metrology

Weak trap: c:\remote\settings202009Sep155709.xml (trap used in momentum interforemeter)

currently working : c:\remote\settings202010Sep115935.xml

tight trap: c:\remote\settings202009Sep160046.xml (tune-out trap)

current: c:\remote\settings202010Sep120216.xml

very tight trap: c:\remote\settings201916Sep182721.xml (qd trap)

current:c:\remote\settings202010Sep122012.xml



2020-09-10 1557 bryce
- setting up going to get bec and then get a PAL running

- got BEC in tight trap c:\remote\settings202010Sep160130.xml
- now to find AL freq
- trap shunt done at 20.28 s
- trap turns off at  20.454000
- got a good al pulse going with 1.31MHz (needs opt)
- c:\remote\settings202010Sep175847.xml
- can fit 20 pulses in with 10 ms samp period
- induce a bit of oscillation by using an impulse in the 
- ramped to even tighter trap with 7v quad 
- c:\remote\settings202010Sep181958.xmls
- ramp done at 20.165500, wait to 20.170000 and start pal
- trap off 20.366000
- opt pal freq
4cyc pulse
frequency(MHz) 	A-num(k)
1.71			1.5
1.61			1.3
1.81			0.98,1.3
1.5				2.9,2.3
1.3				1.3
1.6				2.5,2.5


6 cycle pulse
1.61			6.2,6.2,5.3
1.51			9.9,7.2
1.41			5.3,4.9
1.56			6.6,6.8,
1.59			


i fugure that an oscillation that is 20mm on the detector
changes the zeeman frequency by about 10khz so having a width more than a factor of 10 on that should be good

(1/2*const.mhe*(40e-3)^2)/const.h
=8e3 Hz

for 6 cycles the freq FWHM
0.886/(6*(1/1.6e6))
=0.23e6 HZ


for 10 cycles
141khz

will set this then optimize
1.59	9.8,11,11
1.54	10.5,9.2,10.5
1.64	3.86,5.4
1.61	7.0,6.9
1.55	9.09,5.6


decent AL settings for tight trap
1.550 MHz
1v rms
10 cycles
c:\remote\settings202010Sep185018.xml

now need to cause some more oscillation


2020-09-10 2023 bryce
- managed to get a fitter from tune out working
- took a little data to figure out what the trap freq is 

freq_period=[22.54,10;
             26.31,10.05;
             18.16,9.95;
             13.25,9.90;
             31.44,10.1]
         %26.08,10.05;
- comes ut at 922.4(3) Hz which i think is a lab record =)
- osc has a exp time constant of 0.11 s
- 

- taking some data to make at least a gradient plot for this trap
- shooting for 10 point on that curve

- now taking some data for the amplitude dependence
- done
- tight trap settings with oscilation
- c:\remote\settings202010Sep214907.xml

- i think this should be this trap done
- ill leave it taking stability data overnight
- put this in the 


- Good morning, I think the next thing to do is to get the tune out trap set up with some oscillation in it
  - would like oscillation in both X and Y
  - X can be done with a 0.1ms pulse of the shunt current
  - if this does not excite Y enough we can use the kick coil switch with the little coil inside the trap coil (motion excitation coil)
- then get matalb interfacing with labview to scan out a high resolution aliasing plot


20200911 1021 bryce and jacob
- tight trap for kick and drop
c:\remote\settings202011Sep102156.xml
the delay we want to change is 147
<<<<<<< HEAD
- try to edit this file c:\remote\settings202011Sep105121.xml

Success!
	LabView will now scan over c:\remote\settings202011Sep113201.xml and update the ramp time before kick
	Logfile writing with lines of format
		20200911_120941interfacev5,kick_drop_scan,14,max_T,1.00,T_interval,0.20,this_T,0.20
		which is a bit redundant but the file will only be small
	Need to enable 'Use MATLAB Control' option in Advanced tab

	Need to ensure BEC falls from trap relatively centred
	So adjusting push coil application based on kick-free sequence settings202011Sep121732
		TSO 20.366
		settings202011Sep131751 Centre at [-5,-10] (vs [-7.5,-12])
		Not a great improvement, and still minor clipping unfortunately but let's see what happens if made to oscillate?
		Unfortunately it clips pretty badly. 
		Will go for lunch but beforehand set up a scanning run with the standard trap - see how well this performs
		settings202011Sep133543 centred better with push coils
			centre time is different - difference in trap position account for ~1.5 ms change?
		Obtained 141 shots - just over 2 scans of up to 10ms - before LabView crashed! Darn. 

2020-09-14 0852 kieran

Pal settings for very tight trap:

Freq - 1.550
Ampl - 1.0 Vrms
Offset - 0 Vdc
Phase - 90 degrees

Burst:

delay 0.0 ns
cycles 10
phase = 0 degrees
period 2.306 661 800
source external

settings in labview at the start of day
c:\remote\settings202014Sep085305.xml

Tune-out trap (in progress)
c:\remote\settings202014Sep121715.xml

turned nuller back on

getting a decent pal with all 200 pulses mostly visible


Freq - 1.720
Ampl - 0.1 Vrms
Offset - 0 Vdc
Phase - 90 degrees

Burst:

delay 0.0 ns
cycles 10
phase = 0 degrees
period 2.306 661 800
source external


c:\remote\settings202014Sep144008.xml

2020-09-14 1631 Bryce
- KFT got a good tune out settings running c:\remote\settings202014Sep163153.xml
- will opt BEC a bit and set up a outcoupling period scan
- maanged to get a good atom number improvement to 36k counts
- zeeman slower only allows you to opt for ~9min before it gets to hot
- changed ZS fittings for cleaning
  - nylong tubing will allow for cleaning with ethanol and sodium hydroxide
  - set up tests in beakers with ethanol and degreaser
- made some progress in improving flux
- ZS donger was still in
- source keeps dropping out so fustrating
- lost 2nd mot when adjusting push beam, cant see any signal on "0.6 s mot load" or "4s mot load"
- got push beam back and happy on BEC donger
- still cant see any 2nd mot light
- checked flipper mirrors are in right position
- tried adjusting push beam freq and poinging to see 2nd mot signal
- no luck, heading home 
- leaving source to run at high pressure as sometimes this helps it behave better
- argg =(


2020-09-15 1340 bryce
-se up camera on 2nd mot
- got signal 
- maybe it was becasue the trap supply was off
- got upwards of 158 mv on 0.5s mot load 
- 4s mot load pretty low at 24mv
- source is very unstable
- no BEC i think because the 4s load is too low


