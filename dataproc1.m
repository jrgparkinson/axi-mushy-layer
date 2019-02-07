clear; close all;

r_nums = [25 30 35 40 45 50];

for i = 1:numel(r_nums)
    
    filename = strcat('mushyLayerPrevSteadyState', num2str(r_nums(i)), '.mat');
    %load('heatedWirePrevSteadyState50.mat')
    load(filename)
    %load('testing.mat')
    
    [dr, dz, r_num, z_num] = meshGridProperties(r, z);
    
    Rm=constants('Rm');
    dr=r(1,2)-r(1,1);
    psirr=zeros(r_num);
    psirr(1,:)=nan;
    psirr(end,:)=nan;
    psirr(:,1)=nan;
    psirr(:,end)=nan;
    psir=psirr;
    for J=2:r_num-1
        psirr(:,J)=(psi(:,J+1)-2*psi(:,J)+psi(:,J-1))/(dr^2);
        psir(:,J)=(psi(:,J+1)-psi(:,J-1))/(2*dr);
    end
    figure(1); contourf(r,z,psir)
    figure(1); contourf(r,z,psirr)
    
    thetarr=zeros(r_num);
    thetarr(1,:)=nan;
    thetarr(end,:)=nan;
    thetarr(:,1)=nan;
    thetarr(:,end)=nan;
    thetar=thetarr;
    for J=2:r_num-1
        thetarr(:,J)=(theta(:,J+1)-2*theta(:,J)+theta(:,J-1))/(dr^2);
        thetar(:,J)=(theta(:,J+1)-theta(:,J-1))/(2*dr);
    end
    figure(3); contourf(r,z,thetar)
    figure(3); contourf(r,z,thetarr)
    
    dz=z(2,1)-z(1,1);
    psizz=zeros(z_num);
    psizz(1,:)=nan;
    psizz(end,:)=nan;
    psizz(:,1)=nan;
    psizz(:,end)=nan;
    psiz=psizz;
    for J=2:z_num-1
        psizz(J,:)=(psi(J+1,:)-2*psi(J,:)+psi(J-1,:))/(dz^2);
        psiz(J,:)=(psi(J+1,:)-psi(J-1,:))/(2*dz);
    end
    %figure(2); contourf(r,z,psiz)
    %figure(2); contourf(r,z,psizz)
    
    thetazz=zeros(z_num);
    thetazz(1,:)=nan;
    thetazz(end,:)=nan;
    thetazz(:,1)=nan;
    thetazz(:,end)=nan;
    thetaz=thetazz;
    for J=2:z_num-1
        thetazz(J,:)=(theta(J+1,:)-2*theta(J,:)+theta(J-1,:))/(dz^2);
        thetaz(J,:)=(theta(J+1,:)-theta(J-1,:))/(2*dz);
    end
    %figure(4); contourf(r,z,thetaz)
    %figure(4); contourf(r,z,thetazz)
    
    L=psirr-psir./(r)  +psizz;%+psizz
    %L = psizz;
    R=-Rm*r.*thetar;
    Diff = (L-R)./(abs(R)+1);
    
    %figure(4); mesh(r(1:end-2, :),z(1:end-2, :),Diff(1:end-2, :))
    %figure(4); mesh(r(4:end, :),z(4:end, :),Diff(4:end, :)./L(4:end, :))
    
    max_diffs(i) = max(max(abs(Diff)));
end

plot(1./(r_nums.^2), max_diffs);
xlabel('(\Delta r)^2');
title('Max fractional difference between \psi_{rr}-\psi_r/r + \psi_{zz} and -Rm r \theta_r');
