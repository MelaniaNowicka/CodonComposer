   M S B A M L   `   `   `  ����   CodonComposer'    !CodonOptimizer.Pages.Optimization   \ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml     Jmc;http://schemas.openxmlformats.org/markup-compatibility/2006     @d2http://schemas.microsoft.com/expression/blend/2008     Z VFirstFloor.ModernUI, Version=1.0.6.0, Culture=neutral, PublicKeyToken=2d21ec3cd074c59a0mui&http://firstfloorsoftware.com/ModernUI  -   $��15(�.��� $��0,0,0,41q�  � RowDefinitions  �� $	��120 ��� $	��120 ��� $��* � � ColumnDefinitions �� $
��Auto �.+��� -   $��OriginalCPBscoreTextBox��$
��Left=�$��25��$��154,316,0,0q�$��2�� �� TextWrapping$
 Wrap=�$	��Top=�$��86�� J� RenderTransformOrigin$ 
0.506,0.68&� �� 
IsReadOnly . �� -   $��LoadORFButton��$��Load ORF��$
��Left=�$��31��$��547,50,-3,0q�$	��Top=�$	��178��$��12(�&   System.Windows.Controls.DataGrid  -   $��OptimizationDataGrid��$��2��$	��179��$��20,110,-556,0q�$	��Top=�$
��1270��   CanUserReorderColumns .     CanUserSortColumns .     CanUserResizeRows .     
IsReadOnly . $��1q����	#FF333333$
��Left=�	   SelectionMode$	 Single=��� -   $��OptimizeORFButton��$��Optimize ORF��$
��Left=�$��30��$��20,54,-555,0q�$��2��$	��Top=�$
��1270���� -   
 ��Name$
 UploadRankingButton��$��Upload ranking from file��$
��Left=�$��30��$��
249,50,0,0q�$	��Top=�$	��201��$��12(� � ToolTip _� .���� .+��� $
��Bold"�Upload CC ranking�� 20Please upload a CC ranking generated previously.�� 64Be sure that it is .csv file with the right ranking!�� -   $��CurrentRankingCheckBox��$��Use current CC ranking��$
��Left=�$��26��$��	20,54,0,0q�$	��Top=�$	��152��{�.  $��12(��� -    ��Name$ MaxCPBscoreTextBox��$
��Left=�$��25��$��463,393,0,0q�$��2��$
 Wrap=�$	��Top=�$��87�� . �� -	   $ MinCPBscoreTextBox��$
��Left=�$��25��$��464,435,0,0q�$��2��$
 Wrap=�$	��Top=�$��87�� . $ 
0.494,1.44&��� -
   $��SaveToFASTAButton��$��Save to .FASTA��$
��Left=�$��1117,316,-553,0q�$��2��$	��Top=�$	��173���� $��OR��$
��Left=�$��
191,58,0,0q�$	��Top=�$
��Bold"�$ 0.5,0.05&�$��25��$��20��$��13(��� $	��AND��$
��Left=�$��
486,55,0,0q�$	��Top=�$
��Bold"�$ 0.5,0.05&�$��31��$��20��$��13(��� $$��STEP 1: UPLOAD RANKING AND ORF��$
��Left=�$��20��$	��Top=�$	��284��$
��Bold"�$��0,10,0,0q�$��14(��� $&�� STEP 2: CHOOSE OPTIMIZATION MODE��$
��Left=�$��20��$��0,10,0,0q�$	��Top=�$	��284��$
��Bold"�$��1��$��14(��� $��STEP 3: OPTIMIZE ORF��$
��Left=�$��20��$��0,10,0,0q�$	��Top=�$	��284��$
��Bold"�$��2��$��14(��� $��OR��$
��Left=�$��20��$��
174,61,0,0q�$��1��$	��Top=�$��26��$��13(�$
��Bold"��� $��ORIGINAL CPB:��$
��Left=�$��
20,321,0,0q�$��2��$	��Top=�$
��Bold"�$��14(��� $��LAST MAXIMIZED CPB:��$
��Left=�$��285,399,0,0q�$��2��$	��Top=�$
��Bold"�$	��149��$��14(��� $��LAST MINIMIZED CPB:��$
��Left=�$��285,441,0,0q�$��2��$	��Top=�$
��Bold"�$	��149��$��14(�$ 
-0.997,0.5&��� -   $
 RemoveSelectedButton��$��Remove selected��$
��Left=�$��1117,363,-552,0q�$��2��$	��Top=�$	��173��� -   $��OptimizationProgressBar��$
��Left=�$��6��$��20,89,-554,0q�$��2��$	��Top=�$
��1270��	� -   $��MaximalizeRadioButton��$��Maximize CPB��$
��Left=�$��20��$��	20,61,0,0q�$��1��$	��Top=�$	��133��	� -   $��MinimalizeRadioButton��$��Minimize CPB��$
��Left=�$��20��$��
219,61,0,0q�$��1��$	��Top=�$	��163���� -   $ NcScoreTextBox��$��25��$��924,309,-247,0q�$��2��$
 Wrap=�$	��Top=� . �� $!��Effective number of codons:��$
��Left=�$��606,315,-46,0q�$��2��$	��Top=�$
��Bold"�$	��198��$��14(�$ -1.309,0.763&��� $��BEST MAXIMIZED CPB:��$
��Left=�$��286,319,0,0q�$��2��$	��Top=�$
��Bold"�$	��148��$��14(��� $��BEST MINIMIZED CPB:��$
��Left=�$��286,356,0,0q�$��2��$	��Top=�$
��Bold"�$	��149��$��14(�$ 
-0.997,0.5&��� $-��'Effective number of codons (optimized):��$
��Left=�$��606,352,-117,0q�$��2��$	��Top=�$
��Bold"�$	��272��$��14(�$ -1.309,0.763&��� -   $ NcScoreOptimizedTextBox��$
��Left=�$��25��$��924,352,-250,0q�$��2��$
 Wrap=�$	��Top=�$��86�� . �� -   $ BestMaxCPBscoreTextBox��$
��Left=�$��25��$��464,316,0,0q�$��2��$
 Wrap=�$	��Top=�$��86��$ 
0.506,0.68&� . �� -   $ BestMinCPBscoreTextBox��$
��Left=�$��25��$��464,352,0,0q�$��2��$
 Wrap=�$	��Top=�$��86��$ 
0.506,0.68&� . 