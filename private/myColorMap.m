function map = myColorMap(n,varargin)

if n == 0 || isempty(n)
    map = [];
    return;
end

Red         = [1.0000 0.3137 0.2784];
Yellow      = [1.0000 0.8431 0.0000];
Green       = [0.1804 0.5451 0.3412];
Blue        = [0.3922 0.5843 0.9294];
Purple      = [0.6000 0.1961 0.8000];
DPurple     = [0.192 0.055 0.263];

if nargin == 2 && strcmp(varargin{1},'v1')
    cols = [Red;Yellow;Green;Blue;Purple];
    
    if n > 5
        x = linspace(1,5,n);
    else
        x = 1:n;
    end
    
    map = zeros(n,3);
    for i = 1:3
        map(:,i) = interp1(1:size(cols,1),cols(:,i),x);
    end
elseif nargin == 2 && strcmp(varargin{1},'v1rand')
    cols = [Red;Yellow;Green;Blue;Purple];
    
    if n > 5
        x = linspace(1,5,n);
    else
        x = 1:n;
    end
    
    map = zeros(n,3);
    for i = 1:3
        map(:,i) = interp1(1:size(cols,1),cols(:,i),x);
    end
    
    ind = randperm(n);
    map = map(ind,:);
else
    switch n
        case 1,  map = Red;
        case 2,  map = [Red; Yellow];
        case 3,  map = [Red; Yellow; Green];
        case 4,  map = [Red; Yellow; Green; Blue];
        case 5,  map = [Red; Yellow; Green; Blue; Purple];
        case 6,  map = [Red; Yellow; Green; Blue; Purple; (Purple+DPurple)/2];
        case 7,  map = [Red; Yellow; Green; Blue; (Blue+Purple)/2; Purple; (Purple+DPurple)/2];
        case 8,  map = [Red; Yellow; Green; (Green+Blue)/2; Blue; (Blue+Purple)/2; Purple; (Purple+DPurple)/2];
        case 9,  map = [Red; Yellow; (Yellow+Green)/2; Green; (Green+Blue)/2; Blue; (Blue+Purple)/2; Purple; (Purple+DPurple)/2];
        case 10, map = [Red; (Red+Yellow)/2; Yellow; (Yellow+Green)/2; Green; (Green+Blue)/2; Blue; (Blue+Purple)/2; Purple; (Purple+DPurple)/2];
        otherwise
            m = floor(n/5);
            r = rem(n,5);
            RedYellowC = @(x,y) [linspace(Red(1),Yellow(1), x+y+1)', linspace(Red(2),Yellow(2), x+y+1)', linspace(Red(3),Yellow(3), x+y+1)'];
            YellowGreenC = @(x,y) [linspace(Yellow(1), Green(1), x+y+1)', linspace(Yellow(2), Green(2), x+y+1)', linspace(Yellow(3), Green(3), x+y+1)'];
            GreenBlueC = @(x,y) [linspace(Green(1), Blue(1), x+y+1)', linspace(Green(2), Blue(2), x+y+1)', linspace(Green(3), Blue(3), x+y+1)'];
            BluePurpleC = @(x,y) [linspace(Blue(1), Purple(1), x+y+1)', linspace(Blue(2), Purple(2), x+y+1)', linspace(Blue(3), Purple(3), x+y+1)'];
            PurpleDPurpleC = @(x,y) [linspace(Purple(1), DPurple(1), x+y+1)', linspace(Purple(2), DPurple(2), x+y+1)', linspace(Purple(3), DPurple(3), x+y+1)'];
            
            switch r
                case 0
                    RedYellow = RedYellowC(m,0);
                    YellowGreen = YellowGreenC(m,0);
                    GreenBlue = GreenBlueC(m,0);
                    BluePurple = BluePurpleC(m,0);
                    PurpleDPurple = PurpleDPurpleC(m,0);
                case 1
                    RedYellow = RedYellowC(m,0);
                    YellowGreen = YellowGreenC(m,0);
                    GreenBlue = GreenBlueC(m,0);
                    BluePurple = BluePurpleC(m,0);
                    PurpleDPurple = PurpleDPurpleC(m,1);
                case 2
                    RedYellow = RedYellowC(m,0);
                    YellowGreen = YellowGreenC(m,0);
                    GreenBlue = GreenBlueC(m,0);
                    BluePurple = BluePurpleC(m,1);
                    PurpleDPurple = PurpleDPurpleC(m,1);
                case 3
                    RedYellow = RedYellowC(m,0);
                    YellowGreen = YellowGreenC(m,0);
                    GreenBlue = GreenBlueC(m,1);
                    BluePurple = BluePurpleC(m,1);
                    PurpleDPurple = PurpleDPurpleC(m,1);
                case 4
                    RedYellow = RedYellowC(m,0);
                    YellowGreen = YellowGreenC(m,1);
                    GreenBlue = GreenBlueC(m,1);
                    BluePurple = BluePurpleC(m,1);
                    PurpleDPurple = PurpleDPurpleC(m,1);
            end
            
            RedYellow(end,:) = [];
            YellowGreen(end,:) = [];
            GreenBlue(end,:) = [];
            BluePurple(end,:) = [];
            PurpleDPurple(end,:) = [];
            
            map = [RedYellow;
                YellowGreen;
                GreenBlue;
                BluePurple;
                PurpleDPurple];
    end
end

end