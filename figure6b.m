function hfig=figure6b
    % This software is part of the supplementary material of the publication:
    %
    % Eyherabide HG, Samengo I, When and why noise correlations are important in
    % neural decoding, J Neurosci (2013), 33(45): 17921-17936; doi: 10.1523/JNEUROSCI.0357-13.2013
    %
    % Should you use this code, we kindly request you to cite the aforementioned publication.
    %
    % DESCRIPTION:
    %
    % Calculates the increment in the minimum decoding error \Delta \xi_{NI}^{NIP} attainable 
    % by classical NI decoders for all possible values of P(H,L|S1), and constructs panel Fig 6B.
    %
    % INPUT ARGUMENTS:
    %
    % There are no input arguments
    %
    % OUTPUT ARGUMENT:
    %
    % hfig: Figure handle
    %
    % EXAMPLE:
    %
    % plotf6phlds1;
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
    
    
    
    %% Starts estimations of minimum decoding error for different response probabilities conditioned to stimulus S2
    
    derr=zeros(3,1000);             % This vector will contain the estimations of the minimum decoding error
    
    ps1=5;
    
    for indps=1:3,
        ps1=ps1/10;
        ps=[ps1,1-ps1];
        
        for iphlds1=1:1000,
            
            % Sets response probabilities conditioned to stimulus S2
            pr1r2ds=zeros(2,3,3);
            pr1r2ds(1,3,1)=iphlds1/2000;
            pr1r2ds(2,1,1)=pr1r2ds(1,3,1);
            pr1r2ds(1,1,3)=1-pr1r2ds(1,3,1);
            pr1r2ds(2,3,3)=pr1r2ds(1,1,3);
            
            % Calculates joint stimulus-response probabilities
            for is=1:2,
                psr1r2(is,:,:)=pr1r2ds(is,:,:)*ps(is);
            end
            
            % Estimates the minimum information loss
            derr(indps,iphlds1)=100*decerr('nip',psr1r2,1E-10)/ps1;
        end
    end
    
    
    %% Makes plots
    hfig=figure(); % Figure handle
    options={'markersize',10,'LineWidth',5,'color','blue';  % Options for ps1=0.5
        'markersize',6,'LineWidth',4,'color','red';    % Options for ps1=0.05
        'markersize',6,'LineWidth',2,'color','black'}; % Options for ps1=0.005
    
    % Plots estimations of minimum information loss
    for itype=1:3,
        plot((1:1000)/2000,derr(itype,:),'o',options{itype,:});
        hold on;
    end
    
    % Sets legends, axes and title
    legend('P(S_1)=0.5','P(S_1)=0.05','P(S_1)=0.005',4);
    
    xlim([-.05,0.55]);
    ylim([1,105]);
    xlabel('P(H,L|S_1)');
    ylabel('\Delta \xi_{NI}^{NIP} (%)');
    set(gca,'Yscale','log','YTickLabel',{'1','','100'});
    
    title('Figure 6B - J Neurosci (2013), 33(45): 17921-17936');
end
