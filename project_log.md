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














