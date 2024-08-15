TITLE naG
: Na current for axon. No slow inact.
: M.Migliore Jul. 1997

NEURON {
	SUFFIX naG
	USEION na WRITE ina
	RANGE  gbar,thm1,sigm1,km2,sigm2,kh1,kh2,sigh1,sigh2,thh2
}

PARAMETER {
	gbar = 0.010   	(mho/cm2)	:this should be read in
	eps=1e-7 (mV)
	thm1=1	(mV) :this should be read in
	sigm1=4 (mV)
	km2=0.1
	sigm2=13 (mV)
	kh1=0.012
	kh2=0.2
	sigh1=-20 (mV)
	sigh2=3.5 (mV)
	thh2=1 (mV)	:this should be read in
	:mmin=0.02	
	:hmin=0.5	

	ena=50		(mV)            : must be explicitly def. in hoc
	celsius
	v 		(mV)
}


UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
} 

ASSIGNED {
	ina 		(mA/cm2)
	thegna		(S/cm2) 
	alpham
	betam
	alphah
	betah
	minf	(ms)
	mtau	(ms)
	hinf
	htau
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD euler
        thegna = gbar*m*m*m*h
	ina = thegna * (v - ena)
} 

INITIAL {
	trates(v)
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

PROCEDURE trates(vm) {  
        
    alpham = -((v-thm1-eps)/sigm1)/(exp(-(v-thm1-eps)/sigm1)-1)
    betam = km2*exp(-v/sigm2)
	
	mtau = 1/(alpham+betam)
	minf = alpham/(alpham+betam)

	alphah = kh1/exp(-v/sigh1)
    betah = -kh2*(v-thh2)/(exp(-(v-thh2)/sigh2)-1)
	
	htau =  1/(alphah+betah)
    hinf = alphah/(alphah+betah)                 
	
}

        

        

