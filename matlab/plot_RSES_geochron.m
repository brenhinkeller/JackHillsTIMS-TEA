%% Load dataset
    load rses

    % Eliminate known contamination
    t = rses.Age206Pb238U<2000 & rses.Discordance<30;
    rses.Age206Pb238U(t) = NaN;
    rses.r206Pb_238U(t) = NaN;
    rses.Age207Pb206Pb(t) = NaN;
    rses.Age207Pb235U(t) = NaN;

%% Calculate error ellipse semimajor and semiminor axes

    % Determine sigma level at which to plot
    p = 0.05; % p-value
    sigmaLevel = sqrt(chi2inv(1-p,2)); % The bivariate p=0.05 level
    % sigmaLevel = 2; % Naive two-sigma
    
    % Calculate 1-sigma absolute uncertainties for each isotopic ratio
    % (converting from 2-sigma percent)
    r206Pb_238U_sigma = rses.r206Pb_238U .* rses.r206Pb_238U_2sigma / 2 / 100;
    r207Pb_235U_sigma = rses.r207Pb_235U .* rses.r207Pb_235U_2sigma / 2 / 100;
    
    % Allocate output variables
    a = NaN(size(rses.r206Pb_238U));
    b = a;
    beta = a;
    
% Loop through each analysis
    for i=1:length(rses.r206Pb_238U)
        % Construct covariance matrix from uncertainties and correlation coefficients
        covmat = [r206Pb_238U_sigma(i).^2, r206Pb_238U_sigma(i).*r207Pb_235U_sigma(i).*rses.Corr6875(i);
            r206Pb_238U_sigma(i).*r207Pb_235U_sigma(i).*rses.Corr6875(i), r207Pb_235U_sigma(i).^2];
        
        % For analysis with data
        if ~any(isnan(covmat))
            
            % Calculate eigenvectors and eigenvalues from the covariance matrix.
            % V: matrix of eigenvectors, D: diagonal matrix of eigenvectors
            [V,D] = eig(covmat);
            % Find index of major and minor axes
            [~,major] = max(diag(D));
            [~,minor] = min(diag(D));
            % Larger (major) eigenvector
            v = V(:,major);
            
            % Calculate angle of ellipse from horizontal
            beta(i) = atan(v(2)/v(1));
            % Calculate length of semimajor and semiminor axes for given p-value
            a(i) = sigmaLevel*sqrt(abs(D(major,major)));
            b(i) = sigmaLevel*sqrt(abs(D(minor,minor)));
        end
    end

    
%% Plot dates on concordia diagram, colored by leachate number
    
    figure; hold on;
    c = lines(4); % Colors
    
% Start with dots at the center of each ellipse
    for i=1:4
        t = rses.L==i;
        plot(rses.r207Pb_235U(t),rses.r206Pb_238U(t),'.','MarkerSize',10,'Color',c(i,:));
    end
    legend('L1','L2','L3','R')
    
% Add error ellipses
    for i=1:4
        t = rses.L==i;
        ellipse(b(t),a(t),-beta(t),rses.r207Pb_235U(t),rses.r206Pb_238U(t),c(i,:),100);
    end
    
    xlabel('207Pb / 235U'); ylabel('206Pb / 238U');
    formatfigure;
    
%% One zoom option

    xlim([41.5,60.5])
    xticks(42:2:60)
    ylim([0.808,0.912])
    yticks(0.81:0.01:0.91)
    
%% Concordant analyses only
    figure; hold on;
    c = lines(4); % Colors
    
% Start with dots at the center of each ellipse
    for i=2:4
        t = rses.L==i;
        t = t & rses.Discordance<0.7;
        plot(rses.r207Pb_235U(t),rses.r206Pb_238U(t),'.','MarkerSize',10,'Color',c(i,:));
    end
    legend('L2','L3','R')

% Add error ellipses
    for i=2:4
        t = rses.L==i;
        t = t & rses.Discordance<0.7;
        ellipse(b(t),a(t),-beta(t),rses.r207Pb_235U(t),rses.r206Pb_238U(t),c(i,:),100);
    end
    
    xlabel('207Pb / 235U'); ylabel('206Pb / 238U');
    formatfigure;
    
%% Specific zircon only
    
    figure; hold on;
    c = lines(4); % Colors
    
% Start with dots at the center of each ellipse
    for i=1:4
        t = rses.L==i;
        t = t & contains(rses.zircon,'z6.10');
        plot(rses.r207Pb_235U(t),rses.r206Pb_238U(t),'.','MarkerSize',10,'Color',c(i,:));
    end
    legend('L1','L2','L3','R')

% Add error ellipses
    for i=1:4
        t = rses.L==i;
        t = t & contains(rses.zircon,'z6.10');
        ellipse(b(t),a(t),-beta(t),rses.r207Pb_235U(t),rses.r206Pb_238U(t),c(i,:),100);
    end
    
    xlabel('207Pb / 235U'); ylabel('206Pb / 238U');
    formatfigure;
  
    xlim([51.75,59.75])
    ylim([0.8675,0.9095])

%% Add concordia curve with uncertainty

% Jaffey decay constants
    l235U = log(2)/(7.0381*10^8); % 1/Years
    l235U_sigma = log(2)/(7.0381*10^8)*0.0048/7.0381; % 1/Years
    l238U = log(2)/(4.4683*10^9); % 1/Years
    l238U_sigma = log(2)/(4.4683*10^9)*0.0024/4.4683; % 1/Years
    
    % Uncertainty of 235 decay constant relative to the 238 decay constant
    l235U_sigma_r = l235U .* sqrt((l238U_sigma/l238U).^2 + (l235U_sigma/l235U).^2); % 1/Years
    
% Schoene adjusted 235U decay constant
    l235U_Schoene = 9.8569E-10; % 1/Years
        
% Plot the concordia curve
    xl = xlim; % Note current size of figure
    tlim = log(xl+1)./l235U; % Calculate time range of current window
    t = linspace(tlim(1),tlim(2),1000); % Time vector
    r75t = exp(l235U.*t) - 1; % X axis values
    r68t = exp(l238U.*t) - 1; % Y axis values
    hold on; fill([exp((l235U-l235U_sigma_r*2).*t) - 1, fliplr(exp((l235U+l235U_sigma_r*2).*t) - 1)],[r68t, fliplr(r68t)],'k','FaceAlpha',0.2,'EdgeAlpha',0) % Two-sigma concordia uncertainty
    hold on; plot(r75t,r68t,'k') % Concordia line
    
    r75t_Schoene = exp(l235U_Schoene.*t) - 1; % X axis values
    hold on; plot(r75t_Schoene,r68t,'--k') % Concordia line

    xlim(xl); % Ensure that figure size hasn't changed
    
% Calculate desired range of age markers
    scale=floor(log10(range(tlim))); % Characteristic timescale (order of magnitude)
    trange = round(tlim/10.^scale); % Minimum and maximum time to a round number
    majorstep = 0.25;
    tticks = (trange(1):majorstep:trange(2)).*10^scale; % Time ticks, to a round number
    r75tticks = exp(l235U.*tticks) - 1; % X axis values
    r68tticks = exp(l238U.*tticks) - 1; % Y axis values

% Plot age markers with text labels
    hold on; plot(r75tticks,r68tticks,'.k','MarkerSize',12)
    ticklabels = cellfun(@num2str, num2cell(tticks/10^6),'UniformOutput',false);
    xoffset = range(xlim)/200;
    yoffset = range(ylim)/100;
    test = (r75tticks-xoffset)>(xl(1)+8*xoffset) & r75tticks < xl(2);
    hold on; text(r75tticks(test)-xoffset,r68tticks(test)+yoffset,ticklabels(test),'HorizontalAlignment','right')
  
    fig = gcf;
    warning('off', 'MATLAB:print:FigureTooLargeForPage')
    fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
    set(fig,'renderer','painters')
    

%% Compare SIMS and TIMS ages

% Import the sims dataset from tab-delimited CSV
    sims = importdataset('../data/Compiled_SIMS_Ages.csv','\t');

% For each zircon name in rses.zircon, find the index of a matching name in sims.zircon
    t = findmatches(rses.zircon, sims.zircon);

% Assign non-matches to the last row, which is all NaNs
    t(isnan(t)) = length(sims.zircon);

% Transfer all variables from sims to rses
    for i = 1:length(sims.elements)
        rses.(sims.elements{i}) = sims.(sims.elements{i})(t);
    end

% Plot results
    tmin = 3800;
    tmax = 4250;

    figure; hold on;
    c = lines(4); % Colors
    h = [];
    for i=1:4
        t = rses.L==i & rses.Age207Pb206Pb>tmin;
        t = t & rses.Discordance<0.7;
        hh = herrorbar(rses.Age207Pb206Pb(t),rses.SIMS_Age207Pb206Pb(t),2*rses.Age207Pb206Pb_sigma(t),'.k');
        set(hh,'Color',c(i,:));
        h(i) = errorbar(rses.Age207Pb206Pb(t),rses.SIMS_Age207Pb206Pb(t),2*rses.SIMS_Age207Pb206Pb_sigma(t),'.','MarkerSize',10,'Color',c(i,:));
    end
%     legend(h,{'L1','L2','L3','R'})    
    
    plot([tmin tmax],[tmin tmax],'-k')
    xlim([tmin tmax])
    xlabel('TIMS 6/7 Age')
    ylabel('SIMS 6/7 Age')
    formatfigure
    
%% Compare SIMS and TIMS ages

% Plot results
    tmin = 3900;
    tmax = 4200;

    figure; hold on;
    c = lines(4); % Colors
    h = [];
    for i=4:4
        t = rses.L==i & rses.Age207Pb206Pb>tmin;
        t = t & rses.Discordance<0.7;
        hh = herrorbar(rses.Age207Pb206Pb(t),rses.SIMS_Age207Pb206Pb(t),2*rses.Age207Pb206Pb_sigma(t),'.k');
        set(hh,'Color',c(i,:));
        h(i) = errorbar(rses.Age207Pb206Pb(t),rses.SIMS_Age207Pb206Pb(t),2*rses.SIMS_Age207Pb206Pb_sigma(t),'.','MarkerSize',10,'Color',c(i,:));
        
        d = rses.Age207Pb206Pb(t) - rses.SIMS_Age207Pb206Pb(t);
        e = sqrt(rses.SIMS_Age207Pb206Pb_sigma(t).^2 + rses.Age207Pb206Pb_sigma(t).^2);
        wmean(d,e,'print')
    end
%     legend(h,{'L1','L2','L3','R'})    
    

    
    plot([tmin tmax],[tmin tmax],'-k')
    xlim([tmin tmax])
    xlabel('TIMS 6/7 Age')
    ylabel('SIMS 6/7 Age')
    formatfigure

%% Plot RSES leachate pair U/Pb ratio slopes

    figure;
    xlim([0 60])
    xticks(0:10:60)
    ylim([0 1])
    yticks(0:0.2:1)
    
    xlabel('207Pb / 235U'); ylabel('206Pb / 238U');
    AddConcordiaCurveJaffey;
    formatfigure;
    
    l1s = find(rses.L == 1);
    l2s = find(rses.L == 2);
    l3s = find(rses.L == 3);
    residues = find(rses.L == 4);

    % List of 
    myl1 = l1s(findmatches(rses.fragment, rses.fragment(l1s)));
    myl2 = l2s(findmatches(rses.fragment, rses.fragment(l2s)));
    % myl3 = l3s(findmatches(rses.fragment, rses.fragment(l3s)));
    myresidue = residues(findmatches(rses.fragment, rses.fragment(residues)));
    
    no_l3s = isnan(findmatches(rses.fragment, rses.fragment(l3s)));
    lpenultimate = [l3s; myl2(rses.L==2 & no_l3s)];
    mylpenultimate = lpenultimate(findmatches(rses.fragment, rses.fragment(lpenultimate)));

    % Color scheme
    c = lines(4);
    alpha = 0.25;
    x = 0:60;
    
%     t0 = 3950;
%     x0 = exp(l235U.*t0*1E6) - 1;
%     y0 = exp(l238U.*t0*1E6) - 1;

    % L1 - L2 pairs
    t = rses.L==2;
    x0 = rses.r207Pb_235U(t);
    y0 = rses.r206Pb_238U(t);
    delta68 = rses.r206Pb_238U(t) - rses.r206Pb_238U(myl1(t));
    delta75 = rses.r207Pb_235U(t) - rses.r207Pb_235U(myl1(t));
    slope = delta68 ./ delta75;
    
    for i=1:length(slope)
        hold on; plot(x, slope(i)*(x-x0(i))+y0(i), 'Color',[c(1,:),alpha])
    end
    slope = nanmedian(slope);
    hold on; plot(x, slope*(x-nanmedian(x0))+nanmedian(y0), 'Color',c(1,:), 'LineWidth', 2)

    % L2 - L3 pairs
    t = rses.L==3;
    x0 = rses.r207Pb_235U(t);
    y0 = rses.r206Pb_238U(t);
    delta68 = rses.r206Pb_238U(t) - rses.r206Pb_238U(myl2(t));
    delta75 = rses.r207Pb_235U(t) - rses.r207Pb_235U(myl2(t));
    slope = delta68 ./ delta75;
    
    for i=1:length(slope)
        hold on; plot(x, slope(i)*(x-x0(i))+y0(i), 'Color',[c(2,:),alpha])
    end
    slope = nanmedian(slope);
    hold on; plot(x, slope*(x-nanmedian(x0))+nanmedian(y0), 'Color',c(2,:), 'LineWidth', 2)

    % L2/3 - R pairs
    t = rses.L==4;
    x0 = rses.r207Pb_235U(t);
    y0 = rses.r206Pb_238U(t);
    delta68 = rses.r206Pb_238U(t) - rses.r206Pb_238U(mylpenultimate(t));
    delta75 = rses.r207Pb_235U(t) - rses.r207Pb_235U(mylpenultimate(t));
    slope = delta68 ./ delta75;

    for i=1:length(slope)
        hold on; plot(x, slope(i)*(x-x0(i))+y0(i), 'Color',[c(3,:),alpha])
    end
    slope = nanmedian(slope);
    hold on; plot(x, slope*(x-nanmedian(x0))+nanmedian(y0), 'Color',c(3,:), 'LineWidth', 2)
     
    
    
%% Plot histogram of RSES leachate pair lower intercepts


    l1s = find(rses.L == 1);
    l2s = find(rses.L == 2);
    l3s = find(rses.L == 3);
    residues = find(rses.L == 4);

    % List of 
    myl1 = l1s(findmatches(rses.fragment, rses.fragment(l1s)));
    myl2 = l2s(findmatches(rses.fragment, rses.fragment(l2s)));
    myresidue = residues(findmatches(rses.fragment, rses.fragment(residues)));
    
    no_l3s = isnan(findmatches(rses.fragment, rses.fragment(l3s)));
    lpenultimate = [l3s; myl2(rses.L==2 & no_l3s)];
    mylpenultimate = lpenultimate(findmatches(rses.fragment, rses.fragment(lpenultimate)));

    % Color scheme
    c = lines(4);
    alpha = 0.25;
        
    % Calculate 1-sigma absolute uncertainties for each isotopic ratio
    % (converting from 2-sigma percent)
    r206Pb_238U_sigma = rses.r206Pb_238U .* rses.r206Pb_238U_2sigma / 2 / 100;
    r207Pb_235U_sigma = rses.r207Pb_235U .* rses.r207Pb_235U_2sigma / 2 / 100;
    
    % Allocate output variables
    a = NaN(size(rses.r206Pb_238U));
    b = a;
    beta = a;
    
    LowerIntercept = [];
    
    % L1 - L2 pairs
    t = rses.L==2 ...
        & ~isnan(rses.r207Pb_235U) & ~isnan(rses.r206Pb_238U) ...
        & ~isnan(rses.r207Pb_235U(myl1)) & ~isnan(rses.r206Pb_238U(myl1));
    for n = 1:1000        
        for i=find(t)'
            mu_l2 = [rses.r207Pb_235U(i), rses.r206Pb_238U(i)];
            cov_l2 = r206Pb_238U_sigma(i).*r207Pb_235U_sigma(i).*rses.Corr6875(i);
            covmat_l2 = [r206Pb_238U_sigma(i).^2, cov_l2;
                      cov_l2, r207Pb_235U_sigma(i).^2];
            r = mvnrnd(mu_l2,covmat_l2);
            r75_l2 = r(1);
            r68_l2 = r(2);

            mu_l1 = [rses.r207Pb_235U(myl1(i)), rses.r206Pb_238U(myl1(i))];
            cov_l1 = r206Pb_238U_sigma(myl1(i)).*r207Pb_235U_sigma(myl1(i)).*rses.Corr6875(myl1(i));
            covmat_l1 = [r206Pb_238U_sigma(myl1(i)).^2, cov_l1;
                      cov_l1, r207Pb_235U_sigma(myl1(i)).^2];
            r = mvnrnd(mu_l1,covmat_l1);
            r75_l1 = r(1);
            r68_l1 = r(2);

            delta68 = r68_l2 - r68_l1;
            delta75 = r75_l2 - r75_l1;
            slope = delta68 ./ delta75;

            li = fzero(@(x) slope * (exp(l235U.*x) - 1 - r75_l1) + r68_l1 - exp(l238U.*x) + 1, 0);
            LowerIntercept = [LowerIntercept, li];

        end
    end
    
    figure; hist(LowerIntercept, 100)
    
    
    LowerIntercept = [];
    
    % L2 ? L3 pairs
    t = rses.L==3 ...
        & ~isnan(rses.r207Pb_235U) & ~isnan(rses.r206Pb_238U) ...
        & ~isnan(rses.r207Pb_235U(myl2)) & ~isnan(rses.r206Pb_238U(myl2));
    for n = 1:1000        
        for i=find(t)'
            mu_l3 = [rses.r207Pb_235U(i), rses.r206Pb_238U(i)];
            cov_l3 = r206Pb_238U_sigma(i).*r207Pb_235U_sigma(i).*rses.Corr6875(i);
            covmat_l3 = [r206Pb_238U_sigma(i).^2, cov_l3;
                      cov_l3, r207Pb_235U_sigma(i).^2];
            r = mvnrnd(mu_l3,covmat_l3);
            r75_l3 = r(1);
            r68_l3 = r(2);

            mu_l2 = [rses.r207Pb_235U(myl2(i)), rses.r206Pb_238U(myl2(i))];
            cov_l2 = r206Pb_238U_sigma(myl2(i)).*r207Pb_235U_sigma(myl2(i)).*rses.Corr6875(myl2(i));
            covmat_l2 = [r206Pb_238U_sigma(myl2(i)).^2, cov_l2;
                      cov_l2, r207Pb_235U_sigma(myl2(i)).^2];
            r = mvnrnd(mu_l2,covmat_l2);
            r75_l2 = r(1);
            r68_l2 = r(2);

            delta68 = r68_l3 - r68_l2;
            delta75 = r75_l3 - r75_l2;
            slope = delta68 ./ delta75;

            li = fzero(@(x) slope * (exp(l235U.*x) - 1 - r75_l2) + r68_l2 - exp(l238U.*x) + 1, 1E6);
            LowerIntercept = [LowerIntercept, li];

        end
    end
    
    figure; hist(LowerIntercept, 100)
    
        LowerIntercept = [];
    
    % L2/3 ? R pairs
    t = rses.L==3 ...
        & ~isnan(rses.r207Pb_235U) & ~isnan(rses.r206Pb_238U) ...
        & ~isnan(rses.r207Pb_235U(mylpenultimate)) & ~isnan(rses.r206Pb_238U(mylpenultimate));
    for n = 1:1000        
        for i=find(t)'
            mu_l4 = [rses.r207Pb_235U(i), rses.r206Pb_238U(i)];
            cov_l4 = r206Pb_238U_sigma(i).*r207Pb_235U_sigma(i).*rses.Corr6875(i);
            covmat_l4 = [r206Pb_238U_sigma(i).^2, cov_l4;
                      cov_l4, r207Pb_235U_sigma(i).^2];
            r = mvnrnd(mu_l4,covmat_l4);
            r75_l4 = r(1);
            r68_l4 = r(2);

            mu_lpenultimate = [rses.r207Pb_235U(mylpenultimate(i)), rses.r206Pb_238U(mylpenultimate(i))];
            cov_lpenultimate = r206Pb_238U_sigma(mylpenultimate(i)).*r207Pb_235U_sigma(mylpenultimate(i)).*rses.Corr6875(mylpenultimate(i));
            covmat_lpenultimate = [r206Pb_238U_sigma(mylpenultimate(i)).^2, cov_lpenultimate;
                      cov_lpenultimate, r207Pb_235U_sigma(mylpenultimate(i)).^2];
            r = mvnrnd(mu_lpenultimate,covmat_lpenultimate);
            r75_lpenultimate = r(1);
            r68_lpenultimate = r(2);

            delta68 = r68_l4 - r68_lpenultimate;
            delta75 = r75_l4 - r75_lpenultimate;
            slope = delta68 ./ delta75;

            li = fzero(@(x) slope * (exp(l235U.*x) - 1 - r75_lpenultimate) + r68_lpenultimate - exp(l238U.*x) + 1, 1E6);
            LowerIntercept = [LowerIntercept, li];

        end
    end
    
    figure; hist(LowerIntercept, 100)
    
