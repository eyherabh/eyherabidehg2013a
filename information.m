function info=information(psr1r2)
    % This software is part of the supplementary material of the publication:
    %
    % Eyherabide HG, Samengo I, When and why noise correlations are important in
    % neural decoding, J Neurosci (2013), 33(45): 17921-17936; doi: 10.1523/JNEUROSCI.0357-13.2013
    %
    % Should you use this code, we kindly request you to cite the aforementioned publication.
    %
    % DESCRIPTION:
    %
    % Calculates the Shannon mutual information between the stimulus and the
    % population response.
    %
    % This function is currently limited to populations of 2 neurons.
    %
    % INPUT ARGUMENTS:
    %
    %   - psr1r2: the JOINT stimulus-response probability distribution
    %   (P(R1,R2,S), not P(R1,R2|S)). It must be a matrix of size KxNxM, with
    %   K stimuli, N possible values for responses elicited by neuron R1
    %   and M possible values for responses elicited by neuron R2.
    %
    % EXAMPLE:
    %
    % psr1r2(1,1,2)=.25;
    % psr1r2(1,2,1)=.25;
    % psr1r2(2,2,2)=.25;
    % psr1r2(2,3,3)=.25;
    %
    % information(psr1r2);
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
    
    ps=sum(sum(psr1r2,2),3);
    pr1r2=sum(psr1r2);
    
    info=entropy(ps)+entropy(pr1r2)-entropy(psr1r2);
end

