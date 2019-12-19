classdef HMMdotEM
% This class is used to configure HMMdotEM
  
properties(Constant)
  % Replace this part with whatever your path is.
  % for example:
  % PATH = 'C:\Documents and Settings\User\My Documents\Matlab\HMMdotEM';
  % PATH = '/home/user/HMMdotEM';
  PATH = HMMdotEM.deducehmmdotempath;
  
  %% Directory names
  CONFIGDIR = 'config';
  TESTDIR = 'test'; % Where code testing takes place
  BINDIR = 'bin';
  SRCDIR = 'src';
  
  %% Paths
  SRCPATH = [HMMdotEM.PATH filesep HMMdotEM.SRCDIR];
  BINPATH = [HMMdotEM.PATH filesep HMMdotEM.BINDIR];
  TESTPATH = [HMMdotEM.PATH filesep HMMdotEM.TESTDIR];
  CONFIGPATH = [HMMdotEM.PATH filesep HMMdotEM.CONFIGDIR];
  
  %% A useful enumeration
  SETUPPATHS = 0;
  COMPILEMEXFILES   = 1;
  INCLUDETESTS = 2;
  
  %% Mex Files
  MEXFILES = {'fastalphascaled.c','fastbetascaled.c','fastviterbi.c',...
  'fastsamplepath.c'};
end

methods(Static)
  function startup()
    HMMdotEM.configure(HMMdotEM.SETUPPATHS);
  end
  function install()
    HMMdotEM.configure([HMMdotEM.SETUPPATHS HMMdotEM.COMPILEMEXFILES]);
  end
  function configure(configLevels)
    % configuration of the HMMdotEM package
    if nargin==0
      configLevels = [HMMdotEM.SETUPPATHS HMMdotEM.COMPILEMEXFILES];
    end
    fprintf('Running HMMdotEM configuration script...\n');
    for cL = configLevels
      switch cL
        case HMMdotEM.SETUPPATHS
          fprintf('Adding necessary folders to path...\n');
          HMMdotEM.setuppaths;
        case HMMdotEM.COMPILEMEXFILES
          fprintf('Attempting to compile all predefined mex files...\n');
          HMMdotEM.compilemexfiles;
        case HMMdotEM.INCLUDETESTS
          fprintf('Adding testing scripts to path...\n');
          HMMdotEM.includetests;
        otherwise
          error('HMMdotEM:configure:badInput',...
            'bad input');
      end
    end
  end

  function setuppaths
    fprintf('Adding directory:\n\t%s\n',HMMdotEM.CONFIGPATH);
    addpath(genpath(HMMdotEM.CONFIGPATH));
    fprintf('Adding directory:\n\t%s\n',HMMdotEM.SRCPATH);
    addpath(genpath(HMMdotEM.SRCPATH));
    fprintf('Adding directory:\n\t%s\n',HMMdotEM.BINPATH);
    addpath(genpath(HMMdotEM.BINPATH));
  end
  function compilemexfiles
    mF = HMMdotEM.MEXFILES;
    for fI = 1:length(mF)
      f = mF{fI};
      [pDummy mexFun extDummy] = fileparts(f); %#ok<NASGU>
      if ~(exist(mexFun,'file')==3)
        HMMdotEM.mex(f);
      else
        fprintf('%s already on path\n',mexFun);
      end
    end
  end
  function mex(f)
    try
      mex(which(f),'-outdir',HMMdotEM.BINPATH);
      fprintf('mex compiled %s\n',f);
    catch MexExcept %#ok<NASGU>
      fprintf('Could not compile %s\n',f);
    end
  end
  
  function includetests
    fprintf('Adding directory:\n\t%s\n',HMMdotEM.TESTPATH);
    addpath(genpath(HMMdotEM.TESTPATH));
  end
  
  function p = deducehmmdotempath
    pC = which('HMMdotEM'); % Should not be empty
    assert(~isempty(pC),'the classname of "HMMdotEM" has changed');
    index = strfind(pC,[filesep HMMdotEM.CONFIGDIR filesep 'HMMdotEM']);
    i = index(end);
    p = pC(1:i-1);
%     pS = which('startup');
%     pWD = pwd;
    
  end
  function cleanpath(strToken)
    % Remove all paths that contain strToken using strfind
    p = path;
    regStr = [pathsep '*' pathsep];
    pCell = regexp(p,regStr,'split');
    notHasToken = cellfun(@isempty,strfind(pCell,strToken));
    pCell = pCell(notHasToken);
    pStr = sprintf('%s:',pCell{:});
    path(pStr);
  end
  function svnTok = cleanpathsvn
    % Remove all paths that are under svn directories
    % svnTok = cleanpathsvn
    % svnTok probably will be '.svn' (but may be platform specific)

    svnTok = '.svn';
    HMMdotEM.cleanpath(svnTok);
  end
  
  %% Messaging
  function [ss width] = surroundlines(s,c)
    % Surround the string with two columns made up of string c
    % [ss width] = surroundlines(s,c)
    % ss is the result string,
    % width is the number of characters each of its lines
    % If s is empty the result will be [c ' ' c]
    NEWLINE = char(10); % ASCII '\n'==10
    SPACE = ' ';
    if isempty(s)
      s = SPACE;
    end
    % Add a newline to the end to make code cleaner
    if s(end)~=NEWLINE
      s = [s NEWLINE];
    end
    ind = [0 strfind(s,NEWLINE)];
    lineL = ind(2:end)-ind(1:end-1)-1; % The newline is 1 char long
    maxL = max(lineL);
    cL = length(c);
    width = maxL+2*(cL+1); % A space will be added to both sides
    ss = '';
    numLines = length(ind)-1;
    for i=1:numLines
      ss = [ss c SPACE]; %#ok<AGROW>
      temp = [s((ind(i) + 1):(ind(i+1) - 1))...
        repmat(SPACE,1,maxL-lineL(i))];
      ss = [ss temp SPACE c NEWLINE]; %#ok<AGROW>
    end
    ss = ss(1:end-1);
  end
  function s = marknewlines(s,c)
    % Change all occurences of '\n' in s to '\n c '
    NEWLINE = char(10);
    r = sprintf('\n %s ',c);
    s = strrep(s,NEWLINE,r);
  end  

  function s = printversion
    c = '%%';
    v = HMMdotEM.VERSION;
    m = HMMdotEM.MONTH;
    y = HMMdotEM.YEAR;
    a = HMMdotEM.AUTHOR;
    s = sprintf('This is HMMdotEM version %d.%d.%d\nCopyright %s (%02d/%d)\n',...
      v(1),v(2),v(3),a,m,y);
    licenseLink = ['<a href="', HMMdotEM.LICENSEPATH, '">Link to license</a>.'];
    s = [s sprintf('%s ("%s" %s):\n%s\n%s\n',...
      'Please read the license',HMMdotEM.LICENSEFILE,'file',...
      licenseLink,...
      'Running HMMdotEM.printlicense may also work.')];
    if nargout==0
      message = HMMdotEM.surroundlines(s,c);
      fprintf('%s\n',message);
      clear s;
    end
  end
  function printlicense
    type(HMMdotEM.LICENSEPATH);
  end
end
properties(Constant)
  % PLEASE REFER TO THE LICENSE (license.txt) FOR DETAILS,
  % ABOUT THIS INFORMATION AND YOUR RIGHTS.
  VERSION = [0 6 3]; % Also in ReleaseNotes.txt
  MONTH   = 02;
  YEAR    = 2013; % Also in license.txt
  DATE    = [HMMdotEM.MONTH HMMdotEM.YEAR];
  
  AUTHOR = 'Nikola Karamanov';
  
  LICENSEFILE = 'license.txt';
  LICENSEPATH = [HMMdotEM.PATH filesep HMMdotEM.LICENSEFILE];
end
end
