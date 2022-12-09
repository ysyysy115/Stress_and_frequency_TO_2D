function K=readstiffness
    cd('E:\temp')
    copyfile Job-1_STIF1.mtx Job-1_STIF1.txt
    K0 = readmatrix('Job-1_STIF1.txt');
    for iii=1:size(K0,1)
        K(K0(iii,2),K0(iii,3),K0(iii,1)+1)=K0(iii,4);
    end
