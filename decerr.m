function derr=decerr(errtype,psr1r2,precision)
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
    %   errtype: can take for values
    %       'nil': calculates the minimum information loss attainable by NI decoders.
    %       'nip': calculates the minimum information loss attainable by NI decoders 
    %       which estimation criterion is based on NI posterior probabilities.
    %       'map': calculates the information loss induced by a decoder
    %       constructed using the maximum posterior criterion.
    %
    %   psr1r2: the JOINT stimulus-response probability distribution
    %   (P(R1,R2,S), not P(R1,R2|S)). It must be a matrix of size 2xNxM, with
    %   only 2 stimuli, N possible values for responses elicited by neuron R1
    %   and M possible values for responses elicited by neuron R2.
    %
    %   precision: used only when errtype='nil' or 'nip'. Two responses are considered
    %   different only if their noise-independent likelihoods (in the case of 'nil') 
    %   or their noise-independent posteriors (in the case of 'nip') differ in more than
    %   the value specified in "precision".
    %
    % EXAMPLE:
    %
    % psr1r2(1,1,2)=.25;
    % psr1r2(1,2,1)=.25;
    % psr1r2(2,2,2)=.25;
    % psr1r2(2,3,3)=.25;
    %
    % decerr('nil',psr1r2,1E-3);
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
    if ndims(psr1r2)~=3, error('The first argument should be a three dimensional matrix'); end
    
    % Normalise probability if necessary
    normpsr1r2=sum(psr1r2(:));
    if normpsr1r2~=1, psr1r2=psr1r2/normpsr1r2; end
    
    %% Constructing the noise independent probabilities
    
    % Calculates the marginal probability of each neurons and the stimulus, that
    % is "p(Ri,S)", "i" being the index of the neuron.
    for ir=1:2, psr{ir}=squeeze(sum(psr1r2,4-ir)); end
    
    % Calcualtes the probability of the stimulus
    ps=squeeze(sum(psr{1},2));
    
    % Calculates the probability of the responses
    pr1r2=squeeze(sum(psr1r2));
    
    % Calculates the posterior marginal probability, that is, "p(Ri|S)", "i"
    % being the index of the neuron.
    for ir=1:2, prds{ir}=diag(1./ps)*psr{ir}; end
    
    % Calculates the noise-independent likelihood, that is, "pNI(R|S)", where
    % "R=[R1,R2]".
    for is=1:2, pirds(is,:,:)=prds{1}(is,:)'*prds{2}(is,:); end
    
    % Calculates the noise-independent posterior probability, that is,
    % "pNI(S|R)", where "R=[R1,R2]"
    for is=1:2, pisdr(is,:,:)=pirds(is,:,:)*ps(is); end
    pir1r2=sum(pisdr);
    for is=1:2, pisdr(is,:,:)=pisdr(is,:,:)./pir1r2; end
    
    k=isnan(pisdr);
    pisdr(k)=0;
    
    %% Calculation of the information loss
    switch errtype
		case 'min'
			% Calculates the minimum decoding error attainable by a decoder constructed with knowledge of correlations
			derr=sum(min([psr1r2(1,:),psr1r2(2,:)]'));
			
        case 'nil'
            % Minimum decoding error attainable by an NI decoder
            
            psrnil=zeros(size(psr1r2));
            numr=numel(pr1r2);
            
            for ir=1:numr,
                k= (abs(pirds(1,ir)-pirds(1,:))<precision) & (abs(pirds(2,ir)-pirds(2,:))<precision);
                
                for is=1:2,
                    psrnil(is,ir)=mean(psr1r2(is,k));
                end
            end
            
            derr=sum(min([psrnil(1,:);psrnil(2,:)]));
            
       case 'nip'
            % Minimum decoding error attainable by an NI decoder which estimation stage is based on NI posterior probabilities
            
            psrnip=zeros(size(psr1r2));
            numr=numel(pr1r2);
            
            for ir=1:numr,
                k= (abs(pisdr(1,ir)-pisdr(1,:))<precision) & (abs(pisdr(2,ir)-pisdr(2,:))<precision);
                
                for is=1:2,
                    psrnip(is,ir)=mean(psr1r2(is,k));
                end
            end
            
            derr=sum(min([psrnip(1,:);psrnip(2,:)]));
            
        case 'map',
            % Builds the decoded stimulus using the maximum posterior criterion
            % over the noise-independent posterior distribution.
            
            sdecmap=2-(pisdr(1,:,:)>pisdr(2,:,:));
            
            k=find(squeeze(pisdr(1,:,:)==pisdr(2,:,:)) & pr1r2>0);
            derr=inf;
            if isempty(k),
                derr=derrdec(sdecmap);
            else
                for is=1:2,
                    sdecmap(k)=is;
                    derrnew=derrdec(sdecmap);
                    if derr>derrnew, derr=derrnew; end
                end
            end
    end
    
    function dedec=derrdec(sdec)
        % Calculates the decoding error associated with a specific decoder
        % The code is limited to two stimuli.
        
        dedec=sum(psr1r2(1,sdec==2))+sum(psr1r2(2,sdec==1));
    end
end

