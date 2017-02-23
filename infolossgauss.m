function di=infolossgauss(ditype,ps1,mean1,mean2,rho1,rho2)
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
    %   - ditype: can take for values
    %       'map': calculates the information loss induced by a decoder
    %       constructed using the maximum posterior criterion.
    %       'nil': calculates the minimum information loss attainlable by NI
    %       decoders.
    %       'd'  : estimates the minimum information loss attainable by NI
    %       decoders using the Kullback-Leibler divergence between the real and
    %       the noise-independent posterior probabilities.
    %       'dl' : estimates teh minimum information loss attainable by NI
    %       decoders using the Kullback-Leibler divergence between the real
    %       posterior probabilities and the noise-independent posterior
    %       probabilities of order theta.
    %
    %   - ps1: Probability P(S1) of stimulus S1. It must lie between 0 and 1.
    %
    %   - mean1: Vector with two comonents representing the mean value of
    %   responses to stimulus S1 for neuron R1 (first component) and neuron R2
    %   (second component).
    %
    %   - mean2: Vector with two comonents representing the mean value of
    %   responses to stimulus S2 for neuron R1 (first component) and neuron R2
    %   (second component).
    %
    %   - rho1: Correlation coefficient of population responses to stimulus S1. It
    %   must lie between -1 and 1.
    %
    %   - rho2: Correlation coefficient of population responses to stimulus S2. It
    %   must lie between -1 and 1.
    %
    % EXAMPLE:
    %
    % infolossgauss('nil',0.5,[4,4],[6,6],-0.5,0.5);
    %
    % VERSION CONTROL
    %
    % V1.000 Hugo Gabriel Eyherabide, University of Helsinki (20 Nov 2013)
    %
    % Should you find bugs, please contact either Prof. In�s Samengo (samengo at
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
    
    
    %% Calculation of the information loss
    switch ditype
        case 'map'
            % Builds the decoded stimulus using the maximum posterior criterion
            % over the noise-independent posterior distribution.
            
            % Calculates probabilities pss(S,S^{Dec}) of stimulus S and decoded stimulus S^{Dec}.
            K1=mean2(2)-mean1(2);
            K2=(log(ps1/(1-ps1))+(sum(mean2.^2)-sum(mean1.^2))/2)/K1;
            K1=(mean2(1)-mean1(1))/K1;
            
            pss(2,2)=integral2(@(x,y)probsr(1-ps1,x,y,mean2,c2),-Inf,Inf,@(x)(K2-K1*x),Inf);
            pss(2,1)=(1-ps1)-pss(2,2);
            pss(1,1)=integral2(@(x,y)probsr(ps1,x,y,mean1,c1),-Inf,Inf,-Inf,@(x)(K2-K1*x));
            pss(1,2)=ps1-pss(1,1);
            
            % Calculates stimulus-response mutual information
            info=informationgauss(ps1,mean1,mean2,rho1,rho2);
            
            % Calculates MAP information loss
            di=info-(-sum(pss)*log2(sum(pss)')-sum(pss,2)'*log2(sum(pss,2))+pss(:)'*log2(pss(:)));
        case 'nil'
            % Calculates stimulus-response mutual information
            info=informationgauss(ps1,mean1,mean2,rho1,rho2);
            
            % Calculates stimulus entropy
            ents=-ps1*log2(ps1)-(1-ps1)*log2(1-ps1);
            
            % Calculates stimulus entropy conditioned to the response
            entsdr=integral2(@(x,y)entsdrterm(ps1,x,y,mean1,c1,mean2,c2),-Inf,Inf,-Inf,Inf,'method','iterated');
            
            di=info-ents+entsdr;
            
        case 'd',
            di=integral2(@(x,y)diultb(ps1,x,y,mean1,c1,mean2,c2,1),-Inf,Inf,-Inf,Inf,'Method','iterated');
        case 'dl',
            [x,di]=fminsearch(@(b)integral2(@(x,y)diultb(ps1,x,y,mean1,c1,mean2,c2,b),-Inf,Inf,-Inf,Inf,'Method','iterated'),1);
    end
    
    
    
    function [pis1dr,pis2dr]=probsdrib(p,x,y,m1,m2,theta)
        % Unnormalised noise independent probability with parameter theta (Eq.
        % 9b in Eyherabide & Samengo (2013).
        %
        % NOTE that variances are assumed to be unity
        
        
        % Calculates noise independent probability of stimulus given response
        pis1dr=p*exp(-theta/2*((x-m1(1)).^2+(y-m1(2)).^2));
        pis2dr=(1-p)*exp(-theta/2*((x-m2(1)).^2+(y-m2(2)).^2));
        
        pisum=pis1dr+pis2dr;
        k=pisum>0;
        
        pis1dr(k)=pis1dr(k)./pisum(k);
        pis2dr(k)=pis2dr(k)./pisum(k);
    end
    
    
    function psr=probsr(p,x,y,m,c)
        % Normalised stimulus-response probability
        
        % Subtract mean values
        xr=x-m(1);
        yr=y-m(2);
        
        % Calculates normalised probability
        psr=(p/c.k)*exp(-0.5*(xr.^2*c.inv(1,1)+yr.^2*c.inv(2,2)+2*(xr.*yr)*c.inv(2,1)));
    end
    
    
    function di=diultb(p,x,y,m1,c1,m2,c2,b)
        % Calculated the contribution of response (x,y) to the information loss
        % (Eq. 9a in Eyherabide & Samengo (2013)).
        %
        % NOTE that variances are assumed to be unity.
        
        % Calculates noise independent stimulus probabilities given responses
        [pis1dr,pis2dr]=probsdrib(p,x,y,m1,m2,b);

        % Calculates stimulus-response probabilities
        ps1r=probsr(p,x,y,m1,c1);
        ps2r=probsr(1-p,x,y,m2,c2);
        psum=ps1r+ps2r;
        
        % Calculates contributions to information loss
        di=zeros(size(ps1r));
        k=pis1dr>0 & ps1r>0;
        di(k)=ps1r(k).*log2(ps1r(k)./psum(k)./pis1dr(k));
        
        k=pis2dr>0 & ps2r>0;
        di(k)=di(k)+ps2r(k).*log2(ps2r(k)./psum(k)./pis2dr(k));
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
    
    function entsdr=entsdrterm(p,x,y,m1,c1,m2,c2)
        % Calculates teh contribution to the conditional entropy of the stimulus
        % given the response of responses (x,y).
        
        % Calculates stimulus probability given the response
        [ps1dr,ps2dr]=pnilsdr(p,x,y,m1,c1,m2,c2);
        
        % Calculates response probability
        pr=probsr(p,x,y,m1,c1)+probsr(1-p,x,y,m2,c2);
        
        % Calcualtes conditional entropy
        entsdr=zeros(size(x));
        
        k=ps1dr>0;
        entsdr(k)=ps1dr(k).*log2(ps1dr(k));
        
        k=ps2dr>0;
        entsdr(k)=entsdr(k)+ps2dr(k).*log2(ps2dr(k));
        
        
        entsdr=-pr.*entsdr;
    end
    
end

