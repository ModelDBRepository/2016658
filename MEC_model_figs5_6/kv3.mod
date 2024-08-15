TITLE Kv3

NEURON {
	SUFFIX kv3
	USEION k WRITE ik
	RANGE  gbar, thn1, sign1, kn2, sign2
}

PARAMETER {
	gbar = 0.010   	(mho/cm2)	:this should be read in
	thn1=1 (mV) :this should be read in
	sign1=12 (mV)
	kn2=0.001 
	sign2=-8.5 (mV)

	ek=-90		(mV)            : must be explicitly def. in hoc
	celsius

	v 		(mV)
}


UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
} 

ASSIGNED {
	ik 		(mA/cm2)
	thegk		(S/cm2)
	alphan
	betan
}
 

STATE { n }

BREAKPOINT {
        SOLVE states METHOD euler
        thegk = gbar*n*n*n*n
	ik = thegk * (v - ek)
} 

INITIAL {
	trates(v)  

}

DERIVATIVE states {   
        trates(v)      
        n' = alphan*(1-n) - betan*n
}

PROCEDURE trates(vm) {  
        
	alphan = -(v-thn1)/(exp(-(v-thn1)/sign1)-1)
	betan = kn2/exp(-v/sign2)


}