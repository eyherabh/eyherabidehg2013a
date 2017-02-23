function ent=entropy(p)
    % This software is part of the supplementary material of the publication:
    %
    % Eyherabide HG, Samengo I, When and why noise correlations are important in
    % neural decoding, J Neurosci (2013), 33(45): 17921-17936; doi: 10.1523/JNEUROSCI.0357-13.2013
    %
    % Should you use this code, we kindly request you to cite the aforementioned publication.
    %
    % DESCRIPTION:
    %
    % Calculates the entropy of a random variable x with probability
    % distribution p(x).
    %
    % Notice that p should be normalised, that is sum(p)=1.
    %
    % INPUT ARGUMENTS:
    %
    %   p: a probability distribution 
    %
    % EXAMPLE:
    %
    % psr1r2(1,1,2)=.25;
    % psr1r2(1,2,1)=.25;
    % psr1r2(2,2,2)=.25;
    % psr1r2(2,3,3)=.25;
    %
    % entropy(psr1r2);
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
    
    
    % Chooses only positive values
    p=p(p>0);
    
    if isrow(p), p=p'; end
    
    % Calculates the entropy
    ent=-p'*log2(p);
end
