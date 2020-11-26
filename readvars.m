function varargout = readvars(filename,varargin)
%READVARS Creates variables by reading from a file.
%    Use the READVARS function to create variables by reading
%    column-oriented data from a file. READVARS automatically determines
%    the file format from its extension as described below.
%
%    [Var1, Var2, ...] = READVARS(FILENAME) creates variables by reading
%    from a file. The requested outputs of READVARS are a one-to-one
%    mapping of the columns stored in the file.
%      FILENAME can be one of these:
%      - For local files, FILENAME must be a full path that contains
%        a filename and file extension. FILENAME can also be a relative path
%        to the current directory or to a directory on the MATLAB path.
%        For example, to import a file on the MATLAB path
%           opts = READVARS('patients.xls');
%      - For remote files, FILENAME must be a full path using an internationalized
%        resource identifier (IRI). For example, to import a remote
%        file from Amazon S3 cloud specify the full IRI for the file:
%           opts = READVARS('s3://bucketname/path_to_file/my_table.xls');
%        For more information on accessing remote data, see "Read Remote
%        Data" in the documentation.
%
%    [Var1, Var2, ...] = READVARS(FILENAME,OPTS) creates variables by
%    reading from a file stored at LOCATION using the supplied
%    ImportOptions OPTS. OPTS specifies selected variable names, variable
%    types, and other information regarding the location of the data.
%
%    [Var1, Var2, ...] = READVARS(___, Name,Value) specifies additional
%    parameters using one or more name-value pair arguments.
%
%   Name-Value Pairs for both TEXT and SPREADSHEET files:
% --------------------------------------------------------
%   'FileType'             -   Specify the file as 'text' or 'spreadsheet'.
%
%   'Range'                -   The range to consider when detecting the data. 
%                              Range can be any of the following:
%                                 start-cell   -   A string or character vector containing a a column letter
%                                                  and a row number, or a 2 element numeric vector
%                                                  indicating the starting row and column.
%                                 full-range   -   A start and end cell separated by colon, e.g. 'C2:N15', or 
%                                                  a four element numeric vector containing start row, start column,
%                                                  end row, end column, e.g. [2 3 15 13].
%                                 row-range    -   A string or character vector containing a starting row number
%                                                  and ending row number, separated by a colon.
%                                 column-range -   A string or character vector containing a starting columns
%                    		   				    letter and ending column letter, separated by a colon.
%                                 row number   -   A numeric scalar indicating the first row where data is found.      
%
%   'NumHeaderLines'       -   The number of header lines in the file.
%
%   'TreatAsMissing'       -   Text which is used in a file to represent
%                              missing data, e.g. 'NA'.
%   
%   'ExpectedNumVariables' -   The expected number of variables.   
%
%
%   'TextType'             -   'char' or 'string' the type to use for
%                              text variables.
%
%   'DatetimeType'         -   'datetime', 'text', or 'exceldatenum' the
%                              type to use for dates.
%                              Note: 'exceldatenum' not allowed for text
%
%   Name-Value Pairs for TEXT only :
% -----------------------------------
%   'Delimiter'                  -   The delimiter(s) to use in the file.
%
%   'NumHeaderLines'             -   The number of header lines in the
%                                    file.
%
%   'CommentStyle'               -   The style of comments in the file.
%
%   'LineEnding'                 -   The line ending for the file.
%
%   'DateLocale'                 -   The locale used to interpret month and
%                                    day names in datetime text. LOCALE
%                                    must be a character vector or scalar
%                                    string in the form xx_YY. See the
%                                    documentation for DATETIME for more
%                                    information.
%
%   'Encoding'                   -   The text encoding of the file.
%
%   'DurationType'               -   'duration' or 'text' the type to use
%                                    for duration.
%
%   'Whitespace'                 -   Characters to treat as whitespace.
%
%   'DecimalSeparator'           -   Character used to separate the integer
%                                    part of a number from the decimal part
%                                    of the number.
%
%   'ThousandsSeparator'         -   Character used to separate the
%                                    thousands place digits.
%
%   'TrimNonNumeric'             -   Logical used to specify that prefixes
%                                    and suffixes must be removed leaving 
%                                    only the numeric data. 
%
%   'ConsecutiveDelimitersRule'  -   What to do with consecutive delimiters
%                                    that appear in the file.
%                                 * 'split' - each delimiter separates a
%                                             single field
%                                 * 'join'  - groups of delimiters separate
%                                             single fields
%                                 * 'error' - Ignored for detection
%                                             (treated as 'split'), but the
%                                             resulting read will error.
%
%   'LeadingDelimitersRule'      -   What to do with delimiters at the beginning of a line.
%
%   Name-Value Pairs for SPREADSHEETS only :
% ------------------------------------------
%   'Sheet'    -   The sheet from which to detect the data.
%
%   'UseExcel' -   A logical value that specifies whether or not to read the
%                  spreadsheet file using Microsoft(R) Excel(R) for
%                  Windows(R). Set 'UseExcel' to one of these values:
%                    false -  Does not open an instance of Microsoft Excel 
%                             to read the file. This is the default setting.
%                             This setting may cause the data to be 
%                             written differently for files with 
%                             live updates (e.g. formula evaluation or plugins).
%                    true  -  Opens an instance of Microsoft Excel to read 
%                             the file on a Windows system with Excel installed.
%
%   Parameters which are also accepted with import options. These may have
%   slightly different behavior when used with import options:
%
%   Name-Value Pairs supported with Import Options OPTS:
% -------------------------------------------------------
%
%       Text only parameters:
%         'DateLocale' -   Override the locale used when importing dates.
%         'Encoding'   -   Override the encoding defined in import options.
%
%       Spreadsheet only parameters:
%         'Sheet'      -   Override the sheet value in the import options.
%         'UseExcel'   -   Same behavior as READVARS without import options.
%
%
%   See also READCELL, READMATRIX, READTABLE, READTIMETABLE.

%   Copyright 2018-2019 The MathWorks, Inc.

try
    func = matlab.io.internal.functions.FunctionStore.getFunctionByName('readvars');
    [varargout{1:nargout}] = func.validateAndExecute(filename,varargin{:});
catch ME
    throw(ME);
end

