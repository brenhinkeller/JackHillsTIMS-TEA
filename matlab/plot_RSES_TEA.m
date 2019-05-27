%% Load dataset
    load rses

%% Plot Age versus discordance
    figure; hold on; 
    for i=1:4
        plot(rses.Age206Pb238U(rses.L==i),rses.Discordance(rses.L==i),'.','MarkerSize',15)
    end
    legend('L1','L2','L3','R')
    xlabel('Age206Pb238U'); ylabel('Discordance')
    set(gca,'yscale','log')
    
    
    figure; hold on; 
    for i=1:4
        plot(rses.Discordance(rses.L==i),rses.Age207Pb206Pb(rses.L==i),'.','MarkerSize',15)
    end
    legend('L1','L2','L3','R')
    ylabel('Age207Pb206Pb'); xlabel('Discordance')
    set(gca,'xscale','log')
    ylim([3200 4200])
    
    
%% Plot Th/U versus discordance
    test = rses.teaconc>10^4; % Eliminate analyses that are too dilute for acceptable precision
    figure; hold on; 
    for i=1:4
        plot(rses.Discordance(test&rses.L==i),rses.Th(test&rses.L==i)./rses.U(test&rses.L==i),'.','MarkerSize',15)
    end
    legend('L1','L2','L3','R')
    xlabel('Discordance'); ylabel('TEA Th/U')

    figure; hold on;
    for i=1:4
        plot(rses.Discordance(test&rses.L==i),rses.Th_U(test&rses.L==i),'.','MarkerSize',15)
    end
    legend('L1','L2','L3','R')
    xlabel('Discordance'); ylabel('TIMS Th/U')
    
    
%% All elements versus age (concordant only)

    test = rses.teaconc>10^4 & rses.Discordance<1 & rses.Age207Pb206Pb>3500; % Eliminate analyses that are too dilute for acceptable precision
    c = lines(4);
    for j = 1:length(rses.elements)
        if isnumeric(rses.(rses.elements{j})) &&~contains('Age',rses.elements(j)) &&~contains('sigma',rses.elements(j)) &&~contains('Pbu',rses.elements(j)) &&~contains('Uu',rses.elements(j))
            figure; hold on;
            for i=2:4
                plot(rses.Age207Pb206Pb(test&rses.L==i),rses.(rses.elements{j})(test&rses.L==i),'.','MarkerSize',15,'Color',c(i,:))
            end
            legend('L2','L3','R')
            xlabel('Pb-Pb Age (Ma)'); ylabel(rses.elements{j});
        end
    end
    
%% All elements versus discordance relative to their residue

    residues = find(rses.L == 4);
    myresidue = residues(findmatches(rses.fragment, rses.fragment(residues)));
    test = rses.teaconc > 10^4; % Eliminate analyses that are too dilute for acceptable precision
    
    for j = 1:length(rses.elements)
        if isnumeric(rses.(rses.elements{j})) && ~contains('sigma',rses.elements(j)) && ~contains('Pbu',rses.elements(j)) && ~contains('Uu',rses.elements(j))
            figure; hold on;
            for i=1:4
                plot(rses.Discordance(test&rses.L==i),rses.(rses.elements{j})(test&rses.L==i)./rses.(rses.elements{j})(myresidue(test&rses.L==i)),'.','MarkerSize',15)
            end
            legend('L1','L2','L3','R')
            xlabel('Percent Discordance'); ylabel([rses.elements{j} ' relative to residue']);
            set(gca,'yscale','log')
        end
    end
    
%% Plot a Zr/Hf ratio vs time
    Num = 'Zr';
    Den = 'Hf';
    
    test = rses.teaconc>10^4 & rses.Discordance<1 & rses.Age207Pb206Pb>3000; % Eliminate analyses that are too dilute for acceptable precision
%     LREEI = (rses.Dy./rses.Nd)+(rses.Dy./rses.Sm);
%     test = test & LREEI>30;
    
    figure; hold on;
    c = lines(4);
    for i=2:4
        plot(rses.Age207Pb206Pb(test&rses.L==i),rses.(Num)(test&rses.L==i)./rses.(Den)(test&rses.L==i),'.','MarkerSize',15,'Color',c(i,:));
    end
    legend({'L2','L3','R'})
    xlabel('207Pb-207Pb Age (Ma)'); ylabel([Num ' / ' Den]);
    formatfigure;
    
%% Plot LREE-I vs time

    test = rses.teaconc>10^4 & rses.Age207Pb206Pb>3000; % Eliminate analyses that are too dilute for acceptable precision
    figure; hold on;
    c = lines(4); h=zeros(4,1);
    for i=1:4
        h(i) = plot(rses.Age207Pb206Pb(test&rses.L==i),(rses.Dy(test&rses.L==i)./rses.Nd(test&rses.L==i))+(rses.Dy(test&rses.L==i)./rses.Sm(test&rses.L==i)),'.','MarkerSize',15,'Color',c(i,:));
        herrorbar(rses.Age207Pb206Pb(test&rses.L==i),(rses.Dy(test&rses.L==i)./rses.Nd(test&rses.L==i))+(rses.Dy(test&rses.L==i)./rses.Sm(test&rses.L==i)),2*rses.Age207Pb206Pb_sigma(test&rses.L==i),'.k')
    end
    legend(h,{'L1','L2','L3','R'})
    xlabel('207Pb-207Pb Age (Ma)'); ylabel('LREE-I');
    formatfigure;
    
    
%% Plot LREE-I vs discordance

    test = rses.teaconc>10^4 & rses.Age207Pb206Pb>3000; % Eliminate analyses that are too dilute for acceptable precision
    figure; hold on;
    c = lines(4);
    for i=1:4
        plot(rses.Discordance(test&rses.L==i),(rses.Dy(test&rses.L==i)./rses.Nd(test&rses.L==i))+(rses.Dy(test&rses.L==i)./rses.Sm(test&rses.L==i)),'.','MarkerSize',15,'Color',c(i,:));
    end
    legend({'L1','L2','L3','R'})
    xlabel('Percent Discordance'); ylabel('LREE-I');
    set(gca,'XScale','log');
    formatfigure;
    
%% Plot Zr% vs discordance

    test = rses.teaconc>10^4 & rses.Age207Pb206Pb>3000; % Eliminate analyses that are too dilute for acceptable precision
    figure; hold on;
    c = lines(4);
    for i=1:4
        plot(rses.Discordance(test&rses.L==i),rses.Zr(test&rses.L==i)./496000*100,'.','MarkerSize',15,'Color',c(i,:));
    end
    legend({'L1','L2','L3','R'})
    xlabel('Percent Discordance'); ylabel('Zr (%)');
    set(gca,'XScale','log');
    formatfigure;
    
%% Plot LREE-I vs Zr

    test = rses.teaconc>10^4 & rses.Age207Pb206Pb>3000; % Eliminate analyses that are too dilute for acceptable precision
    figure; hold on;
    c = lines(4);
    for i=1:4
        plot(rses.Zr(test&rses.L==i)./496000*100,(rses.Dy(test&rses.L==i)./rses.Nd(test&rses.L==i))+(rses.Dy(test&rses.L==i)./rses.Sm(test&rses.L==i)),'.','MarkerSize',15,'Color',c(i,:));
    end
    legend({'L1','L2','L3','R'})
    xlabel('Zr (%)'); ylabel('LREE-I');
    formatfigure;
    
    figure; hold on;
    c = lines(4);
    for i=1:4
        plot((rses.Dy(test&rses.L==i)./rses.Nd(test&rses.L==i))+(rses.Dy(test&rses.L==i)./rses.Sm(test&rses.L==i)),rses.Zr(test&rses.L==i)./496000*100,'.','MarkerSize',15,'Color',c(i,:));
    end
    legend({'L1','L2','L3','R'})
    ylabel('Zr (%)'); xlabel('LREE-I');
    formatfigure;


%% Plot discordance versus age relative to residue

    residues = find(rses.L == 4);
    myresidue = residues(findmatches(rses.fragment, rses.fragment(residues)));
    test = rses.Discordance<100; % Eliminate analyses that are too dilute for acceptable precision
    
    figure; hold on;
    c = lines(4); h=zeros(4,1);
    for i=1:4
        t = test & rses.L==i;
        h(i) = plot(rses.Discordance(t),rses.Age207Pb206Pb(t)-rses.Age207Pb206Pb(myresidue(t)),'.','MarkerSize',15,'Color',c(i,:));
%         sigma = sqrt(rses.Age207Pb206Pb_sigma(t).^2+rses.Age207Pb206Pb_sigma(myresidue(t)).^2);
%         errorbar(rses.Discordance(t),rses.Age207Pb206Pb(t)-rses.Age207Pb206Pb(myresidue(t)),2*sigma,'.','Color',c(i,:));
    end
    legend(h,{'L1','L2','L3','R'})
    xlabel('Percent Discordance'); ylabel('206Pb-207Pb age relative to residue');
    set(gca,'XScale','log');
    formatfigure;
    
%% Plot spider diagram of leachages relative to residues

    residues = find(rses.L == 4);
    myresidue = residues(findmatches(rses.fragment, rses.fragment(residues)));

%% All elements
%    elem = {'Sc';'Ti';'Mn';'Fe';'Sr';'Y';'Zr';'Nb';'Ba';'La';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Hf';'Ta';'Pbrad';'Th';'U';};
%    plotelem = {'Sc';'Ti';'Mn';'Fe';'Sr';'Y';'Zr';'Nb';'Ba';'La';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Hf';'Ta';'Pb*';'Th';'U';};

% Selected elements
    elem = {'La';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Hf';'Zr';'Y';'Pbrad';'Th';'U'};
    plotelem = {'La';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Hf';'Zr';'Y';'Pb*';'Th';'U'};
    
       
    x = 1:length(elem);

    L1_R = NaN(size(elem));
    L2_R = NaN(size(elem));
    L3_R = NaN(size(elem));
    
    L1_R_l = NaN(size(elem));
    L2_R_l = NaN(size(elem));
    L3_R_l = NaN(size(elem));
    
    L1_R_u = NaN(size(elem));
    L2_R_u = NaN(size(elem));
    L3_R_u = NaN(size(elem));
    for i=1:length(elem)
        L1_R(i) = nanmedian(rses.(elem{i})(rses.L==1)./rses.(elem{i})(myresidue(rses.L==1)));  
        L2_R(i) = nanmedian(rses.(elem{i})(rses.L==2)./rses.(elem{i})(myresidue(rses.L==2)));
        L3_R(i) = nanmedian(rses.(elem{i})(rses.L==3)./rses.(elem{i})(myresidue(rses.L==3)));
        
        L1_R_l(i) = prctile(rses.(elem{i})(rses.L==1)./rses.(elem{i})(myresidue(rses.L==1)),25);
        L2_R_l(i) = prctile(rses.(elem{i})(rses.L==2)./rses.(elem{i})(myresidue(rses.L==2)),25);
        L3_R_l(i) = prctile(rses.(elem{i})(rses.L==3)./rses.(elem{i})(myresidue(rses.L==3)),25);
        
        L1_R_u(i) = prctile(rses.(elem{i})(rses.L==1)./rses.(elem{i})(myresidue(rses.L==1)),75);
        L2_R_u(i) = prctile(rses.(elem{i})(rses.L==2)./rses.(elem{i})(myresidue(rses.L==2)),75);
        L3_R_u(i) = prctile(rses.(elem{i})(rses.L==3)./rses.(elem{i})(myresidue(rses.L==3)),75);
    end
    
    figure; plot(x,L1_R,'.-','MarkerSize',15)
    hold on; plot(x,L2_R,'.-','MarkerSize',15)
    hold on; plot(x,L3_R,'.-','MarkerSize',15)

    figure; errorbar(x,L1_R,L1_R_l,L1_R_u-L1_R,'.-','MarkerSize',15)
    hold on; errorbar(x,L2_R,L2_R-L2_R_l,L2_R_u-L2_R,'.-','MarkerSize',15)
    hold on; errorbar(x,L3_R,L3_R-L2_R_l,L2_R_u-L3_R,'.-','MarkerSize',15)
    set(gca,'xtick',x,'xticklabel',plotelem)
    xlim([0 length(elem)+1])
    ylabel('Median enrichment over zircon residue')
    legend('L1/Residue','L2/Residue','L3/Residue')
    formatfigure
    

    
    
%% Plot mineral partition coeffs relative to zircon
    load partitioncoeffs.mat
    % load partitioncoeffsOrig.mat

    SiO2=70;
    minerals={'Monazite','Allanite','Xenotime','Zircon','Apatite','Sphene','Melt'}; % Accessory


    p.minerals = unique([p.minerals; 'Pyroxene'; 'Melt']);
    p.Pyroxene.elements = p.Orthopyroxene.elements;
    p.Melt.elements = p.Olivine.elements;
    for i=1:length(p.Pyroxene.elements)
        e = p.Pyroxene.elements{i};
        p.Pyroxene.(e) = nanmean([p.Orthopyroxene.(e),p.Clinopyroxene.(e)],2);
        p.Melt.(e) = zeros(size(p.Olivine.(e)));
    end
    minerals=minerals(logical(cellfun(@(x) sum(ismember(p.minerals, x)), minerals)));

    % All elements
    elem = {'Sc';'Ti';'Mn';'Fe';'Sr';'Y';'Zr';'Nb';'Ba';'La';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Hf';'Ta';'Pb';'Th';'U'}; %TEA suite
    % Elements we're relatively confident in    
    % elem = {'Zr';'Nb';'Y';'La';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Hf';'Pb';'Th';'U'};

    % Replace GERM kds with Sano kds for REEs
    p.Zircon.La = log10(0.00046).*ones(size(p.Zircon.La));
    p.Zircon.Ce = log10(0.36).*ones(size(p.Zircon.La));
    p.Zircon.Pr = log10(0.0172).*ones(size(p.Zircon.La));
    p.Zircon.Nd = log10(0.077).*ones(size(p.Zircon.La));
    p.Zircon.Sm = log10(0.80).*ones(size(p.Zircon.La));
    p.Zircon.Eu = log10(1.22).*ones(size(p.Zircon.La));
    p.Zircon.Gd = log10(8.0).*ones(size(p.Zircon.La));
    p.Zircon.Tb = log10(20.7).*ones(size(p.Zircon.La));
    p.Zircon.Dy = log10(45.9).*ones(size(p.Zircon.La));
    p.Zircon.Ho = log10(80).*ones(size(p.Zircon.La));
    p.Zircon.Er = log10(136).*ones(size(p.Zircon.La));
    p.Zircon.Tm = log10(197).*ones(size(p.Zircon.La));
    p.Zircon.Yb = log10(277).*ones(size(p.Zircon.La));
    p.Zircon.Lu = log10(325).*ones(size(p.Zircon.La));


    figure;
    c=lines(length(minerals));
    trend=NaN(length(minerals),length(elem));
    for j=1:length(minerals)
        for i=1:length(elem) 
            if isfield(p.Zircon,elem{i})
    %             trend(j,i)=p.(minerals{j}).(elem{i})(round(SiO2)-40); % Absolute partition coefficients
                trend(j,i)=p.(minerals{j}).(elem{i})(round(SiO2)-40)-p.Zircon.(elem{i})(round(SiO2)-40); 
            end
        end

        ind=1:length(elem);  
        in=ind(~isnan(trend(j,:)));
        tr=trend(j,~isnan(trend(j,:)));

    %     hold on; plot(in,tr,'.-','Color',c(j,:)) % Interpolate through missing elements
        hold on; plot(ind,trend(j,:),'.-','Color',c(j,:),'MarkerSize',10) % Skip missing elements
    end

    set(gca,'XTick',1:length(elem))
    set(gca,'XTickLabel',elem)

    yt=floor(min(get(gca,'ylim'))):ceil(max(get(gca,'ylim')));
    set(gca,'YTick',yt);
    set(gca,'YTickLabel',10.^yt);
    ylabel('Mineral/Melt Partition Coefficient relative to zircon');
    xlim([1 length(elem)])
    legend(minerals)
    formatfigure;
    
%% Plot absolute zircon kds

    elem = {'La';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Hf';'Th';'U'};

    figure;
    j=4;
    trend = NaN(1,length(elem));

        for i=1:length(elem) 
            if isfield(p.Zircon,elem{i})
                trend(i)=p.(minerals{j}).(elem{i})(round(SiO2)-40); % Absolute partition coefficients
            end
        end

        ind=1:length(elem);  
        in=ind(~isnan(trend));
        tr=trend(~isnan(trend));

        hold on; plot(ind,trend,'.-','Color',c(j,:),'MarkerSize',10) % Skip missing elements


        set(gca,'XTick',1:length(elem))
    set(gca,'XTickLabel',elem)

    yt=floor(min(get(gca,'ylim'))):ceil(max(get(gca,'ylim')));
    set(gca,'YTick',yt);
    set(gca,'YTickLabel',10.^yt);
    ylabel('Zircon/Melt Partition Coefficient');
    xlim([1 length(elem)])
    formatfigure;
    
%% Plot median blank versus leaching step

    [x,y,sigma] = binmedians(rses.L,rses.Pbcpg,0.5,4.5,1,4);
    figure; errorbar(x,y,2*sigma,'.-')
    set(gca,'Xtick',1:4,'XTickLabel',{'L1','L2','L3','R'})
    ylabel('Median Pbc (pg)')
    xlim([0.5 4.5])
    setfigurefontsize(14)
    