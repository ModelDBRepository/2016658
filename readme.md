This is the python and NEURON code associated with the paper:
> Upchurch CM, Knowlton CK, Chamberland S Canavier CC,  Persistent Interruption in Parvalbumin Positive Inhibitory Interneurons: Biophysical and Mathematical Mechanisms

This model entry was contributed by C Canavier. The freely available NEURON
simulation enivronment from [nrn.readthedocs.io](https://nrn.readthedocs.io)
and Python is required for this model.

For help downloading and using NEURON models, see [https://modeldb.science/NEURON_Dwnldguide](https://modeldb.science/NEURON_Dwnldguide)


Bifurcations were analyzed using `Matcont`

The equations for our model in `Matcont` are below

```
area=7916.813487046279*1e-8                                            
Iapp=4.5e-7/area #change to 3.3e-7/area for bursting model                                                          
ek=-90                                                                     
taun=(0.087+11.4/(1+exp((V+14.6)/8.6)))*(0.087+11.4/(1+exp(-(V-1.3)/18.7)))
ninf=1/(1+exp(-(V+12.4)/6.8))                                              
minf=(1/(1+exp(-(V+22)/11.5)))                                             
hinf=1/(1 + exp(-(V+58.3)/-6.7))                                           
tauh = 0.5 + 14 / ( 1 + exp(-(V+60)/-12))
Ina=0.1125*M^3*H*(V-50)                                                    
minfa = (1/(1 + exp(-(V+41.4)/26.6)))^4                                    
mtaua =0.5/(3^(1/10))                                                      
Ia=0.005*A*hslowest*(V-ek)                                                 
Ipas=0.00025*(V+65)                                                        
A'=(minfa-A)/mtaua                                                         
M'=(minf-M)/0.001                                                          
H'=(hinf-H)/tauh                                                           
V'=(Iapp-(Ina+Ipas+Ia+0.225*(V-ek)*(N^2)))*1000                            
N'=(ninf-N)/taun
```

For the model Via et al 2022
```
area=8143.766620952326*1e-8                                                          
Iapp=inject*1e-8/area                                                                
ipas=0.0001689986677404316*(v+77.7944307717461)                                      
ina = 0.28767750461978714*m*m*m*h* (v - 50)                                          
alpham = -((v-(-49.87866107497816)-1e-7)/4)/(exp(-(v-(-49.87866107497816)-1e-7)/4)-1)
betam = 0.1*exp(-v/13)                                                               
mtau = 1/(alpham+betam)                                                              
minf = alpham/(alpham+betam)                                                         
alphah = 0.012/exp(-v/-20)                                                           
betah = -0.2*(v-(-53.326527961625075))/(exp(-(v-(-53.326527961625075))/3.5)-1)       
htau =  1/(alphah+betah)                                                             
hinf = alphah/(alphah+betah)                                                         
ikv1 = 0.0009233607616445254*a*a*a*a * (v - (-90))                                   
alphaa = -(v-51.90844000870827)/(exp(-(v-51.90844000870827)/12)-1)                   
betaa = 0.02/exp(-v/-80)                                                             
ikv3 = 0.011065851407902236*n*n*n*n * (v - (-90))                                    
alphan = -(v-10.179873677546377)/(exp(-(v-10.179873677546377)/12)-1)                 
betan = 0.001/exp(-v/(-8.5))                                                         
m'=(minf-m)/mtau                                                                     
h'=(hinf-h)/htau                                                                     
a' = alphaa*(1-a) - betaa*a                                                          
n' = alphan*(1-n) - betan*n                                                          
v' = (Iapp-(ikv1+ikv3+ina+ipas))*1000
  
```