%% Import TEA data

    % Find all the temp .csv files in the current directory
    file = dir('*Hadean*.csv');
    
    tea = struct;
    
    for i=1:length(file)
	% Import the entire csv file as a cell array
        data = importc(file(i).name,',');
        
	% Find the location of useul data in the imported array
        % Column index of elements (look for parenthesis in name
        elem_i = contains(data(:,1),')');
        % Check which rows have data
%         hasData = ~all(cellfun(@isempty, data(elem_i,:)));
        hasData = true;
        idx = (1:size(data,2))-1;
        
        % Row-number constants come from inspection of thermo element2 .csv
        % output files
        % Index of concentrations
        conc_i = ~mod(idx-11,14);
        % Index of relative concentration uncertainties
        rsd_i = ~mod(idx-13,14);
        % Index of sample names
        name_i = ~mod(idx-1,14); name_i(end)=0;
        
	% Construct vector of desired column names
        elements = data(elem_i,1);
        elements = regexprep(elements,'++','pp');
        temp.elements = fieldname([regexprep(elements,'\(.*\)',''); regexprep(elements,'\(.*\)','_sigma');]);
        
	% Construct matrix of data
        temp.data = str2double([data(elem_i,conc_i&hasData); data(elem_i,rsd_i&hasData)])';
        temp.data(temp.data<0)=NaN; % Negative concentrations are not meaningful
        temp = elementify(temp);
        
        temp.analysis = data(1,name_i)';
        temp.elements = [temp.elements; 'analysis'];
        
        tea = concatenatedatasets(tea,temp);
    end
    
    tea.elements = sort(tea.elements); % Clean up elements
    % Find and remove analysis that are either calibration solutions
    % (denoted 'TEA-') or blanks.
    notsamples = contains(tea.analysis,'TEA') | contains(tea.analysis,'Blank');
    tea = unelementify(tea);
    tea.data = tea.data(~notsamples,:);
    tea = elementify(tea,'k');
    
    
%% Convert raw tea concentration dataset into equivalent concentrations in zircon

    load elementmass
    elements = regexprep(tea.elements,'[^a-zA-Z].*',''); % Find just the chemical element abbreviations from the tea struct
    sigmas = contains(tea.elements,'sigma'); % Find uncertainties so we don't confuse them with data
    if ~isfield(tea,'data'); tea = unelementify(tea); end; % Make tea.data matrix if necessary
    
    % New struct to fill
    rses = struct;
    rses.elements = {};
    
    % Fill element-by-element
    for i=1:length(mass.elements)
        % Look for elemental matches
        sameElement = strcmp(elements,mass.elements{i}) & ~sigmas;
        % If we find matches, then combine them together
        if any(sameElement)
            rses.elements = [rses.elements; mass.elements(i)];
            rses.(mass.elements{i}) = nanmean(cell2mat(tea.data(:,sameElement)),2);
        end
    end 
        
    % Normalize to 496000 ppm Zr in zircon
    rses = unelementify(rses);
    teaconc = nansum(rses.data,2); %Total trace element concentration;
    rses.data = bsxfun(@rdivide, rses.data, teaconc/496000);
    rses = elementify(rses);
   
	rses.teaconc = teaconc;
    rses.elements = [rses.elements; 'teaconc'];
   
    rses.analysis = tea.analysis;
    rses.elements = ['analysis'; rses.elements;];
    
    
%% Combine TEA and U-Pb datasets
    
    % Import TIMS U-Pb age dataset
    tims = importdataset('../Compiled_TIMS_Ages.csv',',');
    
    % Add age variables to RSES struct
    t = findmatches(rses.analysis,tims.analysis);
    for i=1:length(tims.elements)
        if isnumeric(tims.(tims.elements{i}))
            rses.elements = [rses.elements; tims.elements(i)];
            rses.(tims.elements{i}) = NaN(size(rses.analysis));
            rses.(tims.elements{i})(~isnan(t)) = tims.(tims.elements{i})(t(~isnan(t)));
        end
    end
    
    % Calculate U concentration from TIMS Th/U since U is ~quantitatively
    % lost on the columns
    rses.U = rses.Th./rses.Th_U;
    rses.Pbc = rses.U.*rses.Pbcpg./rses.Upg;
    rses.Pbrad = rses.U.*rses.Pbradpg./rses.Upg;
    rses.elements = [rses.elements; 'Pbc';'Pbrad'];
    rses.Pb = rses.Pbc+rses.Pbrad;
    

    
%% Fill in additional variables
    rses.fragment = regexprep(rses.analysis,' L[0-9]*$','');
    rses.fragment = regexprep(rses.fragment,' R$','');
    
    rses.zircon = regexprep(rses.fragment,'[a-z]*$','');

    rses.L = NaN(size(rses.analysis));
    rses.L(my_contains(' L1',rses.analysis)) = 1;
    rses.L(my_contains(' L2',rses.analysis)) = 2;
    rses.L(my_contains(' L3',rses.analysis)) = 3;
    rses.L(my_contains(' R',rses.analysis)) = 4;
    
    rses.F = findmatches(rses.fragment,unique(rses.fragment));
    rses.Z = findmatches(rses.zircon,unique(rses.zircon));

    exportdataset(rses,'../RSES.csv',',')