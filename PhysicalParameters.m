function []=PhysicalParameters()
%You should complete this function according to lecture matrials and
%instructions

format long g
    load('matrices_P_and_E.mat')
    P %given projection matrix
    E %given essential matrix    
    %1.physical parameters from P (interior and exterior orientation)
    %2.solve and display projection center (exterior orientation)
    P1=[P(1,1),P(1,2),P(1,3),;
        P(2,1),P(2,2),P(2,3),;
        P(3,1),P(3,2),P(3,3),];
    P2=[P(1,4),P(2,4),P(3,4)]';
    PC=-inv(P1)*P2;
    'Projection Center=',PC    
    PC1=[PC(1,1),PC(2,1),PC(3,1),1];
    
    %3.Solve interior orientation matrix and display corresponding interior orientaion parameters 
    [Y Q]=rq(P1);
    'interior orientation matrix=',Y
    Cx=Y(1,1);
    Cy=Y(2,2);
    a=Y(1,2);
    X0=Y(1,3);
    Y0=Y(2,3);
    
     %4.display interior orientation 
    '2Dtranslation=',Y   
    'Cx and Cy are scale on x or y axes.'
    'Cx=',Cx
    'Cy=',Cy
    ' non-orthogonality,a=',a
    'Principle point=(',X0 ,',' ,Y0 ,')'
    
    %5.Solve and display rotations (exterior orientation)
    I1=[-1,0,0;
       0,-1,0;
       0,0,-1];
    %U=pinv(I)*inv(S)*P;
    %U(4,4)=1;
    %'U=',U 
    %'Rotationmatrix=',R=[U(1,1),U(1,2),U(1,3);
     %                    U(2,1),U(2,2),U(2,3);
      %                   U(3,1),U(3,2),U(3,3)]
    %'Position=', PJ=[-U(1,4),-U(2,4),-U(3,4)]
    R=det(I1*Q)*I1*Q
    'omega=',omega=asin(R(3,2)),omegadegree=omega*180/pi,'degree'
    'phi=',phi=atan(-R(3,1)/R(3,3)),phidegree=phi*180/pi,'degree'
    'kappa=',kappa=atan(-R(1,2)/R(2,2)),kappadegree=kappa*180/pi,'degree'
    
    %'omega=',omega=U(1,1);,atan2(omega,1)*206265/3600,'degree'
    %'phi=',phi=U(2,2);,atan2(phi,1)*206265/3600,'degree'
    %'kappa=',kappa=U(3,3);,atan2(kappa,1)*206265/3600,'degree'
    
    %6.Physical parameters from essential matrix E (relative orientation)
    'Singular Value Decomposition E',[U,S,V]=svd(E)
    W=[0,-1,0;1,0,0;0,0,1]; 
    ZZ=[1,0,0;0,1,0;0,0,0];
    %7.solve R (all) alternatives and angles in degrees
    'R1=',R1=V*W*U'    
    'R2=',R2=V*W'*U'
    E*inv(R1')
    E*inv(R2')
    
    %8.solve the base vector
    'Base vector1=',B=U*ZZ*W*U'
    'Base vector2=',B2=U*ZZ*W'*U'
    
    %9.normalize the length of b
    'Base vector1=',BV=[B(3,2),-B(3,1),B(2,1)]
    'Base vector2=',BV2=[B2(3,2),-B2(3,1),B2(2,1)]
    '||Base vector1||=',BV(1,1)^2+BV(1,2)^2+BV(1,3)^2
    '||Base vector2||=',BV2(1,1)^2+BV2(1,2)^2+BV2(1,3)^2
    %display all alternatives of base vecor b
end

function deg = rad2deg(rad)
deg = rad / pi * 180;
end