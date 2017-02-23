function hfig=plotf1p1g(example)
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
    % decoders when responses are Gaussian-distributed, for all possible values 
    % of the stimulus probabiltiies, and constructs a panel like Fig 1I.
    %
    % The code is designed to work with examples like those shown in Figures 1C,
    % that is, for only two stimuli and two neurons. Gaussian distributions are
    % assumed to have variances equal to unity.
    %
    % The minimum information loss is calculated using four different
    % estimators: \Delta I_{NI}, \Delta I_{NI}^{LS}, \Delta I_{NI}^{D}, and \Delta
    % I_{NI}^{DL}. Because there are only two stimuli, \Delta I_{NI} and \Delta
    % I_{NI}^{LS} coincide (see aforementinoed manuscript).
    %
    % INPUT ARGUMENTS:
    %
    % example: should be a structure with the following fields
    %   mean1: Vector with two comonents representing the mean value of
    %   responses to stimulus S1 for neuron R1 (first component) and neuron R2
    %   (second component).
    %   mean2: Vector with two comonents representing the mean value of
    %   responses to stimulus S2 for neuron R1 (first component) and neuron R2
    %   (second component).
    %   rho1: Correlation coefficient of population responses to stimulus S1. It
    %   must lie between -1 and 1.
    %   rho2: Correlation coefficient of population responses to stimulus S2. It
    %   must lie between -1 and 1.
    %   title:   (OPTIONAL) The title of the figure.
    %   pos:     (OPTIONAL) The position of the legend in the figure. The values can be
    %   integers from 0 to 4. Please see the description in the matlab help for
    %   the function "legend".
    %
    % EXAMPLE:
    %
    % ex.mean1=[4,4];
    % ex.mean2=[6,6];
    % ex.rho1=-0.9;
    % ex.rho2=0.7;
    % plotf1p1g(ex);
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
    if ~isfield(example,'mean1'), errmessage=[errmessage, 'The field "mean1" corresponding to the mean value of responses to stimulus S1 is missing']; end
    if ~isfield(example,'mean2'), errmessage=[errmessage, 'The field "mean2" corresponding to the mean value of responses to stimulus S2 is missing']; end
    if ~isfield(example,'rho1'), errmessage=[errmessage, 'The field "rho1" corresponding to the correlation coefficient of responses to stimulus S1 is missing']; end
    if ~isfield(example,'rho2'), errmessage=[errmessage, 'The field "rho2" corresponding to the correlation coefficient of responses to stimulus S2 is missing']; end
    if ~isempty(errmessage), error('infolossgauss:errorinput',errmessage); end
    
    errmessage=[];
    if ~isvector(example.mean1) || length(example.mean1)~=2, errmessage=[errmessage, 'mean1 must be a vector with two components\n']; end
    if ~isvector(example.mean2) || length(example.mean2)~=2, errmessage=[errmessage, 'mean2 must be a vector with two components\n']; end
    if abs(example.rho1)>1, errmessage=[errmessage, 'The value of rho1 must lie between -1 and 1\n']; end
    if abs(example.rho2)>1, errmessage=[errmessage, 'The value of rho2 must lie between -1 and 1\n']; end
    if ~isempty(errmessage), error('infolossgauss:errorinput',errmessage); end

    if ~isfield(example,'pos'), example.pos=0; end
    
    %% Starts estimations of minimum information loss for different response probabilities conditioned to stimulus S2

    di=zeros(3,39);             % This vector will contain the estimations of the minimum information loss
    ditype={'map','d','dl'};    % List of types of estimators of the minimum informaiton loss

    indp=0;
    for ps1=.025:.025:.999,
        
        indp=indp+1;
        
        % Calculates stimulus-response mutual information
        info=informationgauss(ps1,example.mean1,example.mean2,example.rho1,example.rho2);
        
        % Estimates the minimum information loss
        for itype=1:3,
            di(itype,indp)=100*infolossgauss(ditype{itype},ps1,example.mean1,example.mean2,example.rho1,example.rho2)/info;
        end
    end
         
    %% Makes plots
    hfig=figure(); % Figure handle
    options={'markersize',10,'LineWidth',5,'color','blue';  % Options for \Delta I_{NI}
             'markersize',6,'LineWidth',4,'color','red';    % Options for \Delta I_{NI}^{D}
             'markersize',6,'LineWidth',2,'color','black'}; % Options for \Delta I_{NI}^{DL}

    % Plots estimations of minimum information loss
    for itype=1:3,
        plot(.025:.025:.999,di(itype,:),'o',options{itype,:});
        hold on;
    end
    
    % Sets legends, axes and title
    legend('\Delta I_{NI}','\Delta I_{NI}^{D}','\Delta I_{NI}^{DL}',example.pos);
    
    xlim([-.05,1.05]);
    ylim([-5,55]);
    xlabel('P(S_1)');
    ylabel('Information loss (%)');
    if isfield(example,'title'), title(example.title); end
end
