function di=infoloss(ditype,psr1r2,precision)
    
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
    % the decoder construction.
    %
    % This function is currently limited to 2 stimuli and populations of 2
    % neurons.
    %
    % INPUT ARGUMENTS:
    %
    %   - ditype: can take for values
    %       'nil': calculates the minimum information loss attainable by NI decoders.
    %       'map': calculates the information loss induced by a decoder
    %       constructed using the maximum posterior criterion.
    %       'd'  : estimates the minimum information loss attainable by NI
    %       decoders using the Kullback-Leibler divergence between the real and
    %       the noise-independent posterior probabilities.
    %       'dl' : estimates teh minimum information loss attainable by NI
    %       decoders using the Kullback-Leibler divergence between the real
    %       posterior probabilities and the noise-independent posterior
    %       probabilities of order theta.
    %
    %   - psr1r2: the JOINT stimulus-response probability distribution
    %   (P(R1,R2,S), not P(R1,R2|S)). It must be a matrix of size 2xNxM, with
    %   only 2 stimuli, N possible values for responses elicited by neuron R1
    %   and M possible values for responses elicited by neuron R2.
    %
    %   - precision: used only when errtype='nil' or 'nip'. Two responses are considered
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
    % infoloss('nil',psr1r2,1E-3);
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
    
    %% Calculation of the information loss
    switch ditype
        case {'nil','map','nip'}
            % Calculates entropies of real responses
            entr1r2=entropy(pr1r2);
            entsr1r2=entropy(psr1r2);
            
            switch ditype
                case 'nil'
                    psrnil=zeros(size(psr1r2));
                    numr=numel(pr1r2);
                    
                    for ir=1:numr,
                        k= (abs(pirds(1,ir)-pirds(1,:))<precision) & (abs(pirds(2,ir)-pirds(2,:))<precision);
                        
                        for is=1:2,
                            psrnil(is,ir)=mean(psr1r2(is,k));
                        end
                    end
                    
                    prnil=sum(psrnil);
                    
                    entrnil=entropy(prnil);
                    entsrnil=entropy(psrnil);
                    
                    di=entr1r2-entsr1r2+entsrnil-entrnil;
                    
               case 'nip'
                    psrnip=zeros(size(psr1r2));
                    numr=numel(pr1r2);
                    
                    for ir=1:numr,
                        k= (abs(pisdr(1,ir)-pisdr(1,:))<precision) & (abs(pisdr(2,ir)-pisdr(2,:))<precision);
                        
                        for is=1:2,
                            psrnip(is,ir)=mean(psr1r2(is,k));
                        end
                    end
                    
                    prnip=sum(psrnip);
                    
                    entrnip=entropy(prnip);
                    entsrnip=entropy(psrnip);
                    
                    di=entr1r2-entsr1r2+entsrnip-entrnip;
                    
                case 'map',
                    % Builds the decoded stimulus using the maximum posterior criterion
                    % over the noise-independent posterior distribution.
                    sdecmap=2-(pisdr(1,:,:)>pisdr(2,:,:));
                    
                    k=find(squeeze(pisdr(1,:,:)==pisdr(2,:,:)) & pr1r2>0);
                    di=inf;
                    if isempty(k),
                        di=didec(sdecmap);
                    else
                        for is=1:2,
                            sdecmap(k)=is;
                            dinew=didec(sdecmap);
                            if di>dinew, di=dinew; end
                        end
                    end
            end
        case {'d','dl'}
            
            % Calculates posterior probabilities
            for is=1:2, psdr(is,:,:)=squeeze(psr1r2(is,:,:))./pr1r2; end
            
            ksr1r2=psr1r2>0;
            
            switch ditype
                case 'd',
                    di=didl(1);
                case 'dl',
                    [x,di]=fminsearch(@(x)didl(x),2);
            end
            
    end
    
    function probdl=pdl(theta)
        % Determines the noise independent probability with parameter beta
        %
        % This function also uses the variables "pirds" and "ps" defined in the
        % main body of the function "infoloss"
        
        
        % Unormalised probability
        probdl=zeros(size(pirds));
        
        for is=1:2, 
            k=find(pirds(is,:,:)>0);
            probdl(is,k)=pirds(is,k).^theta*ps(is); 
        end
        
        
        % Normalization
        probultsum=sum(probdl);
        for is=1:2, probdl(is,:,:)=probdl(is,:,:)./probultsum; end
    end
    
    function deltaidl=didl(theta)
        % Calculates the Kullback-Leibler divergence between the real
        % stimulus-response probabilities and the noise-independent
        % stimulus-response probabilities with parameter beta.
        
        pdlsdr=pdl(theta);
        deltaidl=psr1r2(ksr1r2)'*log2(psdr(ksr1r2)./pdlsdr(ksr1r2));
    end
    
    function deltaidec=didec(sdec)
        % Calculates the information loss associated with a specific decoder
        
        % Builds the decoded stimulus probabilities
        pssd=zeros(2,2);
        for is=1:2, for isd=1:2, pssd(is,isd)=sum(psr1r2(is,sdec==isd)); end; end
        
        % Prepares column vectors for calculation of information
        psd=sum(pssd);
        
        % Calculate entropies of decoder
        entsd=entropy(psd);
        entssd=entropy(pssd);
        
        % Calculates information loss induced by the maximum posterior criterion
        deltaidec=entr1r2-entsr1r2+entssd-entsd;
    end
end
