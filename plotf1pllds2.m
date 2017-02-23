function hfig=plotf1pllds2(example)
    % This software is part of the supplementary material of the publication:
    %
    % Eyherabide HG, Samengo I, When and why noise correlations are important in
    % neural decoding, J Neurosci (2013), 33(45): 17921-17936; doi: 10.1523/JNEUROSCI.0357-13.2013
    %
    % Should you use this code, we kindly request you to cite the aforementioned publication.
    %
    % DESCRIPTION:
    %
    % Estimates the minimum information loss attainable by noise independent
    % decoders for all possible values of the response probabilities associated
    % to stimulus S2, and constructs panels like Fig 1J and Fig 1K.
    %
    % The code is designed to work with examples like those shown in Figures 1A
    % and 1B. Particularly, it works assuming that there are only two stimuli.
    %
    % The minimum information loss is calculated using four different
    % estimators: \Delta I_{NI}, \Delta I_{NI}^{LS}, \Delta I_{NI}^{D}, and \Delta
    % I_{NI}^{DL}. Because there are only two stimuli, \Delta I_{NI} and \Delta
    % I_{NI}^{LS} coincide (see aforementinoed manuscript).
    %
    % INPUT ARGUMENTS:
    %
    % example: should be a structure with the following fields
    %   pr1r2ds: The conditional probability of the population response given
    %   the stimulus.
    %   ps:      The stimulus probability. It should be a vector with two
    %   elements.
    %   title:   (OPTIONAL) The title of the figure.
    %   pos:     (OPTIONAL) The position of the legend in the figure. The values can be
    %   integers from 0 to 4. Please see the description in the matlab help for
    %   the function "legend".
    %
    % EXAMPLE:
    % 
    % ex.pr1r2ds(1,:,:)=[0,1/2;1/2,0];
    % ex.ps=[1/4,3/4];
    % plotf1pllds2(ex);
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
    
    if ~isfield(example,'pr1r2ds'), error('The field "pr1r2ds" corresponding to the response probabilities conditioned to the stimulus is missing'); end
    if ~isfield(example,'ps'), error('The field "ps" corresponding to the stimulus probabilities is missing'); end
    if numel(example.ps)~=2, error('The field "ps" should contain only two elements corresponding to only two stimuli'); end
    if ~isfield(example,'pos'), example.pos=0; end
    
    
    % Extracts probabilities
    pr1r2ds=example.pr1r2ds;
    ps=example.ps;
    

    %% Starts estimations of minimum information loss for different response probabilities conditioned to stimulus S2

    di=zeros(3,99);             % This vector will contain the estimations of the minimum information loss
    ditype={'map','d','dl'};    % List of types of estimators of the minimum informaiton loss
    
    for ipllds2=1:99,
        
        % Sets response probabilities conditioned to stimulus S2
        pr1r2ds(2,:,:)=0;
        pr1r2ds(2,1,1)=ipllds2/100;
        pr1r2ds(2,3,3)=1-ipllds2/100;
        
        % Calculates joint stimulus-response probabilities
        for is=1:2,
            psr1r2(is,:,:)=pr1r2ds(is,:,:)*ps(is);
        end
        
        % Calculates stimulus-response mutual information
        info=information(psr1r2);
        
        % Estimates the minimum information loss
        for itype=1:3,
            di(itype,ipllds2)=100*infoloss(ditype{itype},psr1r2)/info;
        end
    end
    
    
    %% Makes plots
    hfig=figure(); % Figure handle
    options={'markersize',10,'LineWidth',5,'color','blue';  % Options for \Delta I_{NI}
             'markersize',6,'LineWidth',4,'color','red';    % Options for \Delta I_{NI}^{D}
             'markersize',6,'LineWidth',2,'color','black'}; % Options for \Delta I_{NI}^{DL}
    
    % Plots estimations of minimum information loss
    for itype=1:3,
        plot((1:99)/100,di(itype,:),'o',options{itype,:});
        hold on;
    end
    
    % Sets legends, axes and title
    legend('\Delta I_{NI}','\Delta I_{NI}^{D}','\Delta I_{NI}^{DL}',example.pos);
    
    xlim([-.05,1.05]);
    ylim([-5,105]);
    xlabel('P(L,L|S_2)');
    ylabel('Information loss (%)');
    
    if isfield(example,'title'), title(example.title); end
end
