      SUBROUTINE GAMMA(XX,GX,IER)                                       
      IF(XX-34.5)6,6,4                                                  
    4 IER=2                                                             
      GX=1.E38                                                          
      RETURN                                                            
    6 X=XX                                                              
      ERR=1.0E-6                                                        
      IER=0                                                             
      GX=1.0                                                            
      IF(X-2.0)50,50,15                                                 
   10 IF(X-2.0)110,110,15                                               
   15 X=X-1.0                                                           
      GX=GX*X                                                           
      GO TO 10                                                          
   50 IF(X-1.0)60,120,110                                               
C        SEE IF X IS NEAR NEGATIVE INTEGER OR ZERO                      
   60 IF(X-ERR)62,62,80                                                 
   62 K=X                                                               
      Y=FLOAT(K)-X                                                      
      IF(ABS(Y)-ERR)130,130,64                                          
   64 IF(1.0-Y-ERR)130,130,70                                           
C        X NOT NEAR A NEGATIVE INTEGER OR ZERO                          
   70 IF(X-1.0)80,80,110                                                
   80 GX=GX/X                                                           
      X=X+1.0                                                           
      GO TO 70                                                          
  110 Y=X-1.0                                                           
      GY=1.0+Y*(-0.5771017+Y*(+0.9858540+Y*(-0.8764218+Y*(+0.8328212+   
     1Y*(-0.5684729+Y*(+0.2548205+Y*(-0.05149930)))))))                 
      GX=GX*GY                                                          
  120 RETURN                                                            
  130 IER=1                                                             
      RETURN                                                            
      END                                                               
