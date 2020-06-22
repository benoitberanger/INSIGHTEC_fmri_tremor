function print_spm_figure( varargin )


if nargin==0, help(mfilename('fullpath')); return; end


%% Check inputs

% Init
%-----
parser = inputParser;

% Set possible paramters
%-----------------------
parser.addParameter(      'model',     [], @(x) is_model                    (x) )
parser.addParameter(       'anat',     [], @(x) is_anat                     (x) )
parser.addParameter(      'coord',     [], @(x) is_coord                    (x) )
parser.addParameter(   'contrast',     [], @(x) is_notempty_cellstr         (x) )
parser.addParameter( 'correction', 'none', @(x) is_correction               (x) )
parser.addParameter(  'threshold',  0.001, @(x) is_notempty_scalar_positive (x) )
parser.addParameter(     'nvoxel',      0, @(x) is_notempty_scalar_positive (x) )
parser.addParameter(     'outdir',     [], @(x) is_notempty_char            (x) )

% Parse (first step of checks)
%------
parser.parse( varargin{:} )

% Harvest
%-------
model      = parser.Results.model     ;
anat       = parser.Results.anat      ;
coord      = parser.Results.coord     ;
contrast   = parser.Results.contrast  ;
correction = parser.Results.correction;
threshold  = parser.Results.threshold ;
nvoxel     = parser.Results.nvoxel    ;
outdir     = parser.Results.outdir    ;

fprintf('Harvested paramters : \n')
disp(parser.Results)


%% Second step of checks

assert( numel(model)==numel(anat), 'numel(model) numel(anat) must be equal')

% Create outdir if necessary
if ~(exist(outdir,'dir')==7)
    mkdir(outdir)
end


%% Lets roll !

for iModel = 1 : numel(model)
    
    disp(model{iModel})
    
    % Initialiaze the result UI
    matlabbatch{1}.spm.stats.results.spmmat = model(iModel);
    matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
    matlabbatch{1}.spm.stats.results.conspec.threshdesc = correction;
    matlabbatch{1}.spm.stats.results.conspec.thresh = threshold;
    matlabbatch{1}.spm.stats.results.conspec.extent = nvoxel;
    matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{1}.spm.stats.results.units = 1;
    matlabbatch{1}.spm.stats.results.export = cell(1, 0);
    spm('defaults','fmri')
    spm_jobman('run',matlabbatch)
    
    % Get some viaraibles generated by the result UI
    SPM  = evalin('base',  'SPM');
    xSPM = evalin('base', 'xSPM');
    hReg = spm_XYZreg('FindReg',spm_figure('GetWin','Interactive'));
    
    for iCon = 1:numel(contrast)
        
        disp(contrast{iCon})
        
        % Get contrast index
        con_index = find( strcmp(contrast{iCon}, {SPM.xCon.name}) );
        if isempty(con_index)
            continue
        end
        
        % Change contrast
        [hReg,xSPM,SPM] = change_contrast( con_index, xSPM );
        
        % Load anat
        spm_sections(xSPM,hReg,anat{iModel})
        
        for iCoord = 1 : size(coord,1)
            
            disp(coord{iCoord,2})
            
            % Set coordinates
            spm_results_ui('SetCoords',coord{iCoord,1})
            
            % Save PNG
            F = spm_figure('GetWin','Graphics');
            [pathstr1,     ~, ~] = fileparts(model{iModel});
            [       ~, name2, ~] = fileparts(pathstr1);
            con_name = contrast{iCon};
            con_name = strrep(con_name,' ','');
            fname = fullfile(outdir,[ coord{iCoord,2} '___' con_name '___' name2]);
            disp(fname)
            saveas(F,fname,'png')
            
        end % iCoord
        
    end % iCon
    
end % iModel


end % function



%% Set contrast

function [hReg,xSPM,SPM] = change_contrast( con_index, xSPM )
xSPM2.swd   = xSPM.swd;
try xSPM2.units = xSPM.units; end
%         xSPM2.Ic    = getfield(get(obj,'UserData'),'Ic');
xSPM2.Ic    = con_index;
if isempty(xSPM2.Ic) || all(xSPM2.Ic == 0), xSPM2 = rmfield(xSPM2,'Ic'); end
xSPM2.Im    = xSPM.Im;
xSPM2.pm    = xSPM.pm;
xSPM2.Ex    = xSPM.Ex;
xSPM2.title = '';
if ~isempty(xSPM.thresDesc)
    if strcmp(xSPM.STAT,'P')
        % These are soon overwritten by spm_getSPM
        xSPM2.thresDesc = xSPM.thresDesc;
        xSPM2.u = xSPM.u;
        xSPM2.k = xSPM.k;
        % xSPM.STATstr contains Gamma
    else
        td = regexp(xSPM.thresDesc,'p\D?(?<u>[\.\d]+) \((?<thresDesc>\S+)\)','names');
        if isempty(td)
            td = regexp(xSPM.thresDesc,'\w=(?<u>[\.\d]+)','names');
            td.thresDesc = 'none';
        end
        if strcmp(td.thresDesc,'unc.'), td.thresDesc = 'none'; end
        xSPM2.thresDesc = td.thresDesc;
        xSPM2.u     = str2double(td.u);
        xSPM2.k     = xSPM.k;
    end
end
hReg = spm_XYZreg('FindReg',spm_figure('GetWin','Interactive'));
xyz  = spm_XYZreg('GetCoords',hReg);
[hReg,xSPM,SPM] = spm_results_ui('setup',xSPM2);
TabDat = spm_list('List',xSPM,hReg);
spm_XYZreg('SetCoords',xyz,hReg);
assignin('base','hReg',hReg);
assignin('base','xSPM',xSPM);
assignin('base','SPM',SPM);
assignin('base','TabDat',TabDat);
% figure(spm_figure('GetWin','Interactive'));

end % function



%% Validation function of parameters

%--------------------------------------------------------------------------
function result = is_notempty_cell( x )

result = ~isempty(x) && iscell(x);

end % function

%--------------------------------------------------------------------------
function result = is_notempty_cellstr( x )

result = is_notempty_cell(x) && iscellstr(x);

end % function

%--------------------------------------------------------------------------
function result = is_notempty_char( x )

result = ~isempty(x) && ischar(x);

end % function

%--------------------------------------------------------------------------
function result = is_notempty_scalar_positive( x )

result = ~isempty(x) && isscalar(x) && abs(x)==x;

end % function

%--------------------------------------------------------------------------
function result = is_model( model )

result = is_notempty_cellstr(model);

for iModel = 1 : numel(model)
    assert( exist(model{iModel},'file')==2 , 'model file does not exist : %s', model{iModel})
end % iModel

end % function

%--------------------------------------------------------------------------
function result = is_anat( anat )

result = is_notempty_cellstr(anat);

for iAnat = 1 : numel(anat)
    assert( exist(anat{iAnat},'file')==2 , 'anat file does not exist : %s', anat{iAnat})
end % iAnat

end % function

%--------------------------------------------------------------------------
function result = is_coord( x )

result = is_notempty_cell( x );

sz = size(x);

if sz(2)>2
    result = 0;
    return
end

if sz(1)>0
    for i = 1 : sz(1)
        result = result && is_vector          (x{i,1});
        result = result && is_notempty_cellstr(x(i,2));
    end
end

end % function

%--------------------------------------------------------------------------
function result = is_correction( correction )

result = is_notempty_char(correction);

assert( any(strcmp(correction, {'FWE','none'})), 'correction can be ''FWE'' or ''none'' ' )


end % function

%--------------------------------------------------------------------------
function result = is_vector( x )

result = ~isempty(x) && isnumeric(x) && isvector(x) && numel(x)==3;

end % function