
        % ===================== update M1 M2 ========================
        
        M1=diag(1./max(sum(abs(U).^2,2).^(1/2),1e-10));
        
        M2=diag(1./max(sum(abs(V).^2,2).^(1/2),1e-10));
        
        % ===================== update sigma ========================
        
        diff = (X-U*V').^2;
        sigma =sqrt(1/(2*mFea)*sum(sum(diff)));
        
        % ===================== update H ========================
        
        H = diag(exp(-1*sum(diff,2)/(2*sigma^2)));
        
        % ===================== update U ========================
        
        XV = H*X*V;   
        VV = V'*V;  
        UVV = H*U*VV;
        if 1
            UM1 = gamma*M1*U;
            UVV = UVV+UM1;
        end
        
        U = U.*(XV./max(UVV,1e-10));    
    
        %===================== update V ========================

        XU = 2*X'*H*U;     
        UU = 2*U'*H*U;  
        VUU = V*UU; 
        Y1=(abs(Y)+Y)/2;
        Y2=(abs(Y)-Y)/2;
        if 1
            XU=XU+2*alpha*W*Z+u*Z+Y2;
            VM = 2*gamma*M2*V;
            VUU =VUU+2*alpha*D*Z+VM+u*V+Y1;
        end
        
        V =V.*(XU./max(VUU,1e-10));
        
        %======================update Z ========================
        
        [r,c]=size(Z);
        for i=1:r
          for j=1:c
            Z(i,j)=max(Z(i,j),0);
          end
        end      
        DCol = full(sum(W,2));
        D = spdiags(DCol,[0],nSmp,nSmp);
        
        L = D - W;
        C=V+(Y/u)-(alpha/u)*L'*V;
        Z=max(C,0);

        %======================update W(S) ========================
        
        P=R+(alpha/2*beta)*V*Z';
        Q=(P+P')/2;
        
        [qsz1,qsz2]=size(Q);
        N=nSmp;
        
       if mod(nIter,2)==1  
        tmp1=(N+ones(1,qsz1)*Q*ones(qsz2,1))/(N*N)*ones(qsz1,1)*ones(1,qsz2);
        tmp2=1/N*Q*ones(qsz2,1)*ones(1,qsz2);
        tmp3=1/N*ones(qsz1,1)*ones(1,qsz1)*Q;
        W=Q+tmp1-tmp2-tmp3;
       end
         
       if mod(nIter,2)==0
        W=max(P,0);
        W=W-diag(diag(W));
       end
        
        %======================update Y ========================
        
        Y=Y+u*(V-Z);
        
        %======================update max ======================
        
        u=max(p*u,1e4);



    





        

        