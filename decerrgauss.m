function decerr=decerrgauss(errtype,ps1,mean1,mean2,rho1,rho2)
    % This software is part of the supplementary material of the publication:
    %
    % Eyherabide HG, Samengo I, When and why noise correlations are important in
    % neural decoding, J Neurosci (2013), 33(45): 17921-17936; doi: 10.1523/JNEUROSCI.0357-13.2013
    %
    % Should you use this code, we kindly request you to cite the aforementioned publication.
    %
    % DESCRIPTION:
    %
    % Estimates the information loss induced by ignoring noise correlations in
    % the decoder construction when responses are Gaussian-distributed.
    %
    % This function is currently limited to 2 stimuli and populations of 2
    % neurons. In addition, variances are assumed to be unity.
    %
    % INPUT ARGUMENTS:
    %
    %   -errtype: can take for values
	%		'min': calculates the minimum decoding error attainable by any decoder.
	%	    'map': calculates the decoding error induced by a NI decoder
    %       constructed using the maximum posterior criterion.
    %       'nil': calculates the minimum decoding error attainlable by NI
    %       decoders.
    %
    %   -ps1: Probability P(S1) of stimulus S1. It must lie between 0 and 1.
    %   mean1: Vector with two comonents representing the mean value of
    %   responses to stimulus S1 for neuron R1 (first component) and neuron R2
    %   (second component).
    %
    %   -mean1: Vector with two comonents representing the mean value of
    %   responses to stimulus S2 for neuron R1 (first component) and neuron R2
    %   (second component).
    %
    %   -rho1: Correlation coefficient of population responses to stimulus S1. It
    %   must lie between -1 and 1.
    %
    %   -rho2: Correlation coefficient of population responses to stimulus S2. It
    %   must lie between -1 and 1.
    %
    % EXAMPLE:
    %
    % decerrgauss('nil',0.5,[4,4],[6,6],-0.5,0.5);
    %
    % VERSION CONTROL
    %
    % V1.000 Hugo Gabriel Eyherabide, University of Helsinki (20 Nov 2013)
    %
    % Should you find bugs, please contact either Prof. Inés Samengo (samengo at
    % cab.cnea.gov.ar) or Hugo Gabriel Eyherabide (hugo.eyherabide at helsinki.fi)
    %
    % LICENSE
    %
    % Copyright 2013 Hugo Gabriel Eyherabide
    %
    % This program is free software: you can redistribute it and/or modify it under
    % the terms of the GNU General Public License as published by the Free Software
    % Foundation, either version 3 of the License, or (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful, but WITHOUT ANY
    % WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    % PARTICULAR PURPOSE. See the GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License along with this
    % program. If not, see http://www.gnu.org/licenses/.
    
    
    
    %% Checking for consistency of the data
    errmessage=[];
    if ps1<0 || ps1>1, errmessage=[errmessage, 'The value of ps1 (first argument) must lie between 0 and 1 \n']; end
    if ~isvector(mean1) || length(mean1)~=2, errmessage=[errmessage, 'mean1 (second argument) must be a vector with two components\n']; end
    if ~isvector(mean2) || length(mean2)~=2, errmessage=[errmessage, 'mean2 (third argument) must be a vector with two components\n']; end
    if abs(rho1)>1, errmessage=[errmessage, 'The value of rho1 (fourth argument) must lie between -1 and 1\n']; end
    if abs(rho2)>1, errmessage=[errmessage, 'The value of rho2 (fifth argument) must lie between -1 and 1\n']; end
    if ~isempty(errmessage), error('infolossgauss:errorinput',errmessage); end
    
    
    %% Transform covariance data to speed up calculations
    % c1 and c2: covariance data with correlations
    % c1i and c2i: covariance data without correlations
    
    [c1,c1i]=transformcov(rho1);
    [c2,c2i]=transformcov(rho2);
    
    function [ck,cki]=transformcov(rhok)
        covk=[1,rhok;rhok,1];
        ck.inv=inv(covk);
        ck.k=2*pi*det(covk)^.5;
        cki.k=2*pi;
    end
    
    
    %% Calculation of the decoding error
    switch errtype
		case 'min'
			% Calculates the minimum decoding error attainable by a decoder constructed with knowledge of correlations
			decerr=integral2(@(x,y)proberrmin(ps1,x,y,mean1,c1,mean2,c2),-Inf,Inf,-Inf,Inf);

        case 'map'
            % Calculates the decision boundary using the NI posterior distribution
            K1=mean2(2)-mean1(2);
            K2=(log(ps1/(1-ps1))+(sum(mean2.^2)-sum(mean1.^2))/2)/K1;
            K1=(mean2(1)-mean1(1))/K1;
            
            % Calculates decoding error of MAP NI decoder
            decerr=1-integral2(@(x,y)probsr(1-ps1,x,y,mean2,c2),-Inf,Inf,@(x)(K2-K1*x),Inf)-integral2(@(x,y)probsr(ps1,x,y,mean1,c1),-Inf,Inf,-Inf,@(x)(K2-K1*x));

        case 'nil'
            
            % Calculates the minimum decoding error achievable by a NI decoder
            decerr=integral2(@(x,y)proberrnil(ps1,x,y,mean1,c1,mean2,c2),-Inf,Inf,-Inf,Inf,'method','iterated');
    end
    
    
    
    
    
    function psr=probsr(p,x,y,m,c)
        % Normalised stimulus-response probability
        
        % Subtract mean values
        xr=x-m(1);
        yr=y-m(2);
        
        % Calculates normalised probability
        psr=(p/c.k)*exp(-0.5*(xr.^2*c.inv(1,1)+yr.^2*c.inv(2,2)+2*(xr.*yr)*c.inv(2,1)));
    end
    
	function perr=proberrmin(p,x,y,m1,c1,m2,c2)
		% Error probability for responses (x,y)
		
		ps1r=probsr(p,x,y,m1,c1);
		ps2r=probsr(1-p,x,y,m2,c2);
		
		perr=zeros(size(x));
		perr(:)=min([ps1r(:),ps2r(:)]');
	end


    function [pnils1dr,pnils2dr]=pnilsdr(p,x,y,m1,c1,m2,c2)
        % Calculates the probability P(S|R^{NIL}) of stimuli S given the representation R^{NIL}
        % of the population response.
        %
        % NOTE that variances are assumed to be unity.
        
        % Calculates responses that are merged
        r1=[x(:)-m1(1),y(:)-m1(2)]; % Original responses
        muvec=m2-m1;
        muvec=muvec/norm(muvec);
        if isrow(muvec), muvec=muvec'; end
        r2=2*(r1*muvec)*muvec'-r1;  % Symmetric responses with respect to the line passing through the mean values
        
        % Calculates conditional probabilities
        pnils1dr=zeros(size(x));
        pnils1dr(:,:)=probsr(p,r1(:,1),r1(:,2),[0,0],c1)+probsr(p,r2(:,1),r2(:,2),[0,0],c1);
        
        pnils2dr=zeros(size(x));
        pnils2dr(:,:)=probsr(1-p,r1(:,1),r1(:,2),m2-m1,c2)+probsr(1-p,r2(:,1),r2(:,2),m2-m1,c2);
        
        pnilr=pnils1dr+pnils2dr;
        k=pnilr>0;
        
        pnils1dr(k)=pnils1dr(k)./pnilr(k);
        pnils2dr(k)=pnils2dr(k)./pnilr(k);
    end
    
    function perr=proberrnil(p,x,y,m1,c1,m2,c2)

        % Error probability for responses (x,y) after the NI assumption.
        
        % Calculates stimulus probability given the response
        [ps1dr,ps2dr]=pnilsdr(p,x,y,m1,c1,m2,c2);
        
        % Calculates response probability
        pr=probsr(p,x,y,m1,c1)+probsr(1-p,x,y,m2,c2);
        
        perr=zeros(size(x));
		perr(:)=min([ps1dr(:),ps2dr(:)]').*pr;
    end
    
end

