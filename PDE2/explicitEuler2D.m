function sol= explicitEuler2D(Nx,Ny,dT,sol,M)
 %Calls the matrix of the system with Euler implemented
sol=M*sol;  %Solution after dt
end