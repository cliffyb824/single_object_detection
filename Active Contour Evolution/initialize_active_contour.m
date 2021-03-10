function [beta0,Eevol] = initialize_active_contour(Sinit,Sprior)

I=Sinit.I; % Log-transformed image
T=Sinit.T;
scalefac=Sinit.scalefac;
cmax=Sinit.cmax;
nrow=Sinit.nrow;

switch Sinit.type
    
    case 'manual' % Initialize as hand segmentations via imfreehand
        
        [~,beta0]=hand_segment(exp(I),T,scalefac,cmax);
        Eevol=[];
        
    case 'prior' % Initialize as shape prior mean with average length and centroid
        
        beta0c=Sprior.muL*Sprior.Sbetamean.betaStd;
        beta0=beta0c+repmat(Sprior.muc,T,1);
        Eevol=[];
        
    case 'coarse-refined'
        
        beta0=find_best_rotation_for_init(Sprior,Sinit);
        Eevol=[];
        
    case 'refined'
        
        delta=Sinit.delta;
        maxit=Sinit.maxit;
        t=Sinit.t;
        
        % Initialize accelerated gradient descent
        
        beta0=find_best_rotation_for_init(Sprior,Sinit);
        % Or don't find best rotation
%         beta0c=Sprior.muL*Sprior.Sbetamean.betaStd;
%         beta0=beta0c+repmat(Sprior.muc,T,1);

        Sbeta0=curve_properties(beta0,t);
        Sinit.lambda([2 4])=0; % Hard coded lambda index
        Eevol=zeros(1,maxit);
        curves{1}=Sbeta0.beta;
        i=1;
        [Eevol(i),Sprior,gradE]=Etotal(Sbeta0,Sinit,Sprior,false,[]);
        lam=0;
        b=compute_similarity_basis(Sbeta0);
        nbasis=length(b);
        beta0y_prev=Sbeta0.beta;

        % keyboard;

        % Active contour accelerated gradient method (AGM) loop 
        while i<maxit %&& error(i)>tol 

            % Update algorithm parameters for AGM
            lam=(1+sqrt(1+4*lam^2))/2;
            gam=(1-lam)/(lam+1);

            % Plot the active contour evolution on top of the image
            if Sinit.toggleplot
                plot_curves_on_image(curves,exp(I),nrow,scalefac,1,2,cmax,'r');
                title(['Iter: ' num2str(i)])
            end
%             frame = getframe(gcf);
%             writeVideo(v,frame);

            % Display the iteration count and current energy in the command window
            if Sinit.display_iter
                disp([num2str(i) ': ' num2str(Eevol(i))]);
            end

            % Project gradient to similarity basis
            gradEproj=zeros(T,2);
            for j=1:nbasis
                gradEproj=gradEproj+inner_prod_q(gradE',b{j}',t)*b{j};
            end

            % Compute intermediate sequence item via gradient descent
            beta0y_cur=Sbeta0.beta-delta*gradEproj;

            % Compute AGM update using intermediate sequence
            beta0=(1-gam)*beta0y_cur+gam*beta0y_prev;
            Sbeta0=curve_properties(beta0,t);
            
            % Update energy, curves, and basis
            i=i+1;    
            [Eevol(i),Sprior,gradE]=Etotal(Sbeta0,Sinit,Sprior,false,[]);
            curves{1}=Sbeta0.beta;
            b=compute_similarity_basis(Sbeta0);
            beta0y_prev=beta0y_cur;
           
        end

        Eevol=Eevol(1:i);
        beta0=Sbeta0.beta;

        % Plot the final solution on top of the image
        if Sinit.toggleplot
            plot_curves_on_image(curves,exp(I),nrow,scalefac,1,2,cmax,'r');
            title(['Iter: ' num2str(i)])
        end

        % Display the iteration count and final energy in the command window
        if Sinit.display_iter
            disp([num2str(i) ': ' num2str(Eevol(i))]);
        end

end




