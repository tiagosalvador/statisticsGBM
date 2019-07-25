function latextable(X,varargin)

% PROGRAM OUTLINE:
% 1 - INITILIZE PROGRAM
% 2 - CONVERT INPUT INTO A SINGLE COLUMN OF '&' DELIMINATED STRINGS
% 3 - ADD VERTICAL HEADINGS TO STRINGS
% 4 - BUILD HEADER ROW(S)
% 5 - COMBINE ROWS AND ADD HORIZONTAL LINES
% 6 - BUILD COLUMN INSTRUCTIONS
% 7 - OUTPUT DATA TO FILE
% 8 - APPLY PREFERENCES
% SUBFUNCTION: subfunction_options
% SUBFUNCTION: convertdata
% APPPYPREF applies the last used directory and property settings
%__________________________________________________________________________

% 1 - INITILIZE PROGRAM
    a = subfunction_options(varargin);

% 2 - CONVERT INPUT INTO A SINGLE COLUMN OF '&' DELIMINATED STRINGS
    r = size(X,1); c = size(X,2);
    A = convertdata(X,a.format);
   
% 3 - ADD VERTICAL HEADINGS TO STRINGS
    Hprfx = '';
    if ~isempty(a.vert) && size(a.vert,1) == r;
        for i = 1:r
            A{i,1} = [sprintf('%s & ',a.vert{i,:}),A{i}];
        end
        Hprfx(1:size(a.vert,2)) = '&';
    end

% 4 - BUILD HEADER ROW(S)
    H = {};
    if ~isempty(a.horiz)&& size(a.horiz,2) == c ; 
        for i = 1:size(a.horiz,1)
            H{i,1} = sprintf('%s & ',a.horiz{i,:}); 
            H{i} = [Hprfx,' ',regexprep(H{i},'& ','\\\',c)];   
        end
    end
    
% 5 - COMBINE ROWS AND ADD HORIZONTAL LINES
    ROWS = [H;A]; h = a.hline; nr = length(ROWS);
    for i = 1:length(h);
        if h(i) == 0;
            ROWS = ['\hline'; ROWS(1:nr)]; 
            h = h + 1; nr = nr + 1; 
        elseif isnan(h(i));
            ROWS{nr+1} = '\hline';
            nr = nr + 1;
        else
            ROWS = [ROWS(1:h(i));'\hline';ROWS(h(i)+1:nr)]; 
            h = h + 1; nr = nr + 1;
        end
    end 

% 6 - BUILD COLUMN INSTRUCTIONS
    col(1:c+size(a.vert,2)) = 'c';
    v = a.vline;
    for i = 1:length(v);
        if v(i) == 0; col = ['|',col];
        elseif isnan(v(i)); col = [col,'|'];
        else
            col = regexprep(col,'c','c|',v(i));
        end
    end

% 7 - OUTPUT DATA TO FILE
    fid = fopen(a.name, 'w'); 
    fprintf(fid,'%s\n',['\begin{tabular}{',col,'}']);
    for i = 1:length(ROWS); fprintf(fid,'%s\n',ROWS{i}); end
    fprintf(fid,'%s\n','\end{tabular}');
    fclose(fid);

% 8 - APPLY PREFERENCES
    applypref(a);

    
%--------------------------------------------------------------------------
% SUBFUNCTION: subfunction_options
function [a] = subfunction_options(in)
% SUBFUNCTION_OPTIONS seperates user inputed property modifiers

% 1 - SET THE DEFAULT SETTINGS 
    % 1.1 - Standard settings
        a.horiz  = {};          a.vert  = {};
        a.hline  = [0,NaN];     a.vline = [];
        a.format = '%3.1f';
        a.name  = '';
        a.save = ''; a.load = '';
        
    % 1.2 - Gather/set user defaults
        if ~ispref('latextable_settings','default');
            setpref('latextable_settings','default',a);
        else
            a = getpref('latextable_settings','default');
        end
    
% 2 - SEPERATE THE DATA FROM OPTIONS
list = fieldnames(a); k = 1; n = length(in);
while k < n
    % 2.1 - Get modifier tag, associated data, and initilize for next loop
        opt = in{k}; value = in{k+1}; match = '';  k = k + 2;

    % 2.2 - Compare modifier tag with available options
        if ischar(opt); 
            match = strmatch(lower(opt),lower(list),'exact');
            if ~isempty(match); a.(list{match}) = value; end
        end
        
    % 2.3 - Load preferences (overwrites existing)   
        if strcmpi(opt,'load') && ispref('latextable_save',value);
            a = getpref('latextable',value);
        end

    % 2.4 - Produce an error message if modifier is not found
        if isempty(match);
            mes = ['The property modifier, ',opt,', was not recoignized.'];
            disp(mes);
        end
end

% 3 - DETERMINE OUTPUT LOCATION
    % 3.1 - Locate the last used directory
    if isempty(a.name);
        if ispref('latextable','lastdir'); 
            loc = getpref('latextable','lastdir');
        else
            loc = cd;
        end
        
    % 3.2 - Prompt the user for the path     
        [nm,pth] = uiputfile('*.tex','Save as...',loc);
        if isnumeric(nm); return; end
        a.name = [pth,nm];
    end
    
% 4 - CORRECT VERTICAL AND HORIZONTAL ARRAYS WITH A SINGLETON DIMENSION
    if size(a.horiz,2) == 1; a.horiz = a.horiz'; end
    if size(a.vert,1) == 1; a.vert = a.vert'; end

%--------------------------------------------------------------------------
function A = convertdata(X,form)
% CONVERTDATA changes numeric and cell arrays into strings for each row
% X must be a numeric array or a cell array, with a string of numeric value
% in each entry.

% Define the size of the array
    r = size(X,1); c = size(X,2);
    
% For completely numeric arrays
    if isnumeric(X);
        for i = 1:r;
            j = 1;
            if sum(strfind(form{j},'e'))
                   A{i,1} = sprintf(['\\num{',form{j},'} & '],X(i,j));
                   A{i,1} = strrep(A{i,1}, 'e+0','e+');
                   A{i,1} = strrep(A{i,1}, 'e-0','e-');
            else
                   A{i,1} = sprintf([form{j},' & '],X(i,j));
            end
            for j = 2:c
               if sum(strfind(form{j},'e'))
                   A{i,1} = horzcat([A{i,1},sprintf(['\\num{',form{j},'} & '],X(i,j))]); 
               else
                   A{i,1} = horzcat([A{i,1},sprintf([form{j},' & '],X(i,j))]);
               end    
            end
            A{i} = regexprep(A{i},'& ','\\\',c);
        end
        
% For cell arrays        
    elseif iscell(X);
        for i = 1:r; 
            A{i,1} = num2str(X{i,1},form);
            for j = 2:c;
                A{i,1} = [A{i,1},' &',num2str(X{i,j},form)];   
            end
            A{i,1} = [A{i,1},'\\'];
        end 
    end
%--------------------------------------------------------------------------
function applypref(a)
% APPPYPREF applies the last used directory and property settings

    % 7.1 - Store the last used directory
        pth = fileparts(a.name);
        if ~ispref('latextable','lastdir');
            addpref('latextable','lastdir',pth);
        else
            setpref('latextable','lastdir',pth);
        end
        
    % 7.2 - Save the options
        if ~isempty(a.save);
            apply = a.save; a.save = [];
            if ~ispref('latextable_settings',apply);
                addpref('latextable_settings',apply,a);
            else
                qst = ['The settings, ',apply,' already exist, do you ',...
                'want to overwrite?'];
                qans = questdlg(qst,'Overwrite?','Yes');
                if strcmpi(qans,'yes'); 
                    setpref('latextable_settings',apply,a); 
                end
            end
        end
