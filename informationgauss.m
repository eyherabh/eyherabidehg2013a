function info=informationgauss(ps1,mean1,mean2,rho1,rho2)
    % This software is part of the supplementary material of the publication:
    %
    % Eyherabide HG, Samengo I, When and why noise correlations are important in
    % neural decoding, J Neurosci (2013), 33(45): 17921-17936; doi: 10.1523/JNEUROSCI.0357-13.2013
    %
    % Should you use this code, we kindly request you to cite the aforementioned publication.
    %
    % DESCRIPTION:
    %
    % Calculates the mutual information between stimulus and response when
    % responses are Gaussian-distributed.
    %
    % This function is currently limited to 2 stimuli and populations of 2
    % neurons. In addition, variances are assumed to be unity.
    %
    % INPUT ARGUMENTS:
    %
    %   - ps1: Probability P(S1) of stimulus S1. It must lie between 0 and 1.
    %
    %   - mean1: Vector with two comonents representing the mean value of
    %   responses to stimulus S1 for neuron R1 (first component) and neuron R2
    %   (second component).
    %
    %   - mean1: Vector with two comonents representing the mean value of
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
    % informationgauss(0.5,[4,4],[6,6],-0.5,0.5);
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
    
    c1=transformcov(rho1);
    c2=transformcov(rho2);
    
    function ck=transformcov(rhok)
        covk=[1,rhok;rhok,1];
        ck.inv=inv(covk);
        ck.k=2*pi*det(covk)^.5;
    end
    
    
    % Calcualtes stimulus entropy
    ents=-ps1*log2(ps1)-(1-ps1)*log2(1-ps1);

    % Calcualtes response entropy
    entr=integral2(@(x,y)entrterm(ps1,x,y,mean1,c1,mean2,c2),-Inf,Inf,-Inf,Inf,'method','iterated');

    % Calcualtes joint stimulus-response entropy
    entsr=integral2(@(x,y)(entsrterm(ps1,x,y,mean1,c1)+entsrterm(1-ps1,x,y,mean2,c2)),-Inf,Inf,-Inf,Inf,'method','iterated');

    % Calculates stimulus-response mutual information
    info=ents+entr-entsr;

    function psr=probsr(p,x,y,m,c)
        % Normalised stimulus-response probability
        
        % Subtract mean values
        xr=x-m(1);
        yr=y-m(2);
        
        % Calculates normalised probability
        psr=(p/c.k)*exp(-0.5*(xr.^2*c.inv(1,1)+yr.^2*c.inv(2,2)+2*(xr.*yr)*c.inv(2,1)));
    end
    
    
    function entr=entrterm(p,x,y,m1,c1,m2,c2)
        % Calculates the contribution to the response entropy of the response
        % (x,y).
        
        % Calcultes the response probability
        pr=probsr(1-p,x,y,m2,c2)+probsr(p,x,y,m1,c1);
        
        % Calculates the response entropy
        k=find(pr>0);
        entr=zeros(size(x));
        entr(k)=-pr(k).*log2(pr(k));
    end
    
    
    function entsr=entsrterm(p,x,y,m,c)
        % Calculates the contribution to the joint stimulus-response entropy of
        % response (x,y) and stimulus with probability p
        
        % Calculates stimulus-response probability
        psrent=probsr(p,x,y,m,c);
        
        % Calculates stimulus-response entropy
        k=find(psrent>0);
        entsr=zeros(size(x));
        entsr(k)=-psrent(k).*log2(psrent(k));
    end

end

