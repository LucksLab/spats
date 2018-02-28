
#import <math.h>

#import "H2Private.h"
#import "UIClient.h"


@interface Profiles : NSObject
@property(strong)    NSArray * treated;
@property(strong)    NSArray * untreated;
@property(strong)    NSArray * beta;
@property(strong)    NSArray * theta;
@property(strong)    NSArray * rho;
@property(assign)    CGFloat c;
-(void)compute;
@end


@interface Site : NSObject
@property(assign) NSInteger L;
@property(assign) NSInteger site;
@property(assign) NSInteger treated;
@property(assign) NSInteger untreated;
@property(assign) CGFloat beta;
@property(assign) CGFloat theta;
@property(assign) CGFloat rho;
-(void)reset;
-(BOOL)valid;
-(Site *)clone;
@end


@interface SelectionPart : NSObject
-(id)initWithSite:(Site *)site;
@property(readonly) Site * site;
@property(readonly) id<H2View> view;
@end


@interface SpatsMatrix : H2ViewImpl_View < H2ViewLayout, H2ViewDrawer, H2MouseTracking >
{
    UIClient * m_client;
    NSInteger m_viewId;
    NSMutableDictionary * m_profiles;
    NSInteger m_minLength;
    NSInteger m_maxLength;
    CGSize m_siteSize;
    CGRect m_matrixFrame;
    Site * m_curSite;

    id<H2View> m_rowTracker;
    id<H2View> m_colTracker;
    id<H2View> m_siteInfo;
    NSMutableArray * m_siteLabels;
    NSMutableArray * m_xLabels;
    NSMutableArray * m_yLabels;

    NSString * m_plotType;
    NSMutableDictionary * m_selectedData;
    CGFloat m_max;

    BOOL m_plotTypeRow;
    NSMutableArray * m_selection;
    BOOL m_showNull;
}

@end


@implementation SpatsMatrix

-(id)init
{
    if ((self = [super init])) {
        self.drawer = self;
        self.layout = self;

        m_curSite = [[Site alloc] init];
        m_siteSize = CGSizeMake(5, 5);
        m_matrixFrame = CGRectMake(0, 0, 800, 800);
        m_plotTypeRow = YES;
        m_showNull = NO;
        [self addClickHandler:[H2EventHandler handlerWithTarget:self selector:@selector(handleClick:)]];

        m_rowTracker = [H2ViewImpl container];
        m_rowTracker.backgroundColor = [H2Color clearColor];
        m_rowTracker.borderColor = [H2Color whiteColor];
        m_rowTracker.borderWidth = 1;
        m_rowTracker.hidden = YES;
        [self addSubview:m_rowTracker];

        m_colTracker = [H2ViewImpl container];
        m_colTracker.backgroundColor = [H2Color clearColor];
        m_colTracker.borderColor = [H2Color whiteColor];
        m_colTracker.borderWidth = 1;
        m_colTracker.hidden = YES;
        [self addSubview:m_colTracker];

        m_siteInfo = [H2ViewImpl container];
        m_siteInfo.backgroundColor = [H2Color darkGrayColor];
        m_siteInfo.cornerRadius = 6;
        m_siteInfo.hidden = YES;
        [self addSubview:m_siteInfo];

        m_siteLabels = [[NSMutableArray alloc] init];
        for (NSInteger i = 0; i < 7; ++i) {
            id<H2Label> l = [H2ViewImpl label];
            l.textColor = [H2Color whiteColor];
            l.fontSize = 14;
            [m_siteInfo addSubview:l];
            [m_siteLabels addObject:l];
        }

    }
    return self;
}

-(void)matrix_plot:(NSDictionary *)params
{
    m_max = [params[@"max"] floatValue];
    m_selectedData = [[NSMutableDictionary alloc] init];
    m_plotType = params[@"plot"];
    SEL prof_sel = NSSelectorFromString(m_plotType);
    for (id key in m_profiles.allKeys) {
        Profiles * p = m_profiles[key];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Warc-performSelector-leaks"
        m_selectedData[key] = [p performSelector:prof_sel];
#pragma clang diagnostic pop
    }
    [self layoutSubviews];
    [self setNeedsDisplay];
}

-(void)show_null
{
    m_showNull = YES;
    [self setNeedsDisplay];
}

-(void)hide_null
{
    m_showNull = NO;
    [self setNeedsDisplay];
}

-(void)row
{
    m_plotTypeRow = YES;
    [self resetSelection];
}

-(void)column
{
    m_plotTypeRow = NO;
    [self resetSelection];
}

-(void)updateWithModel:(NSDictionary *)model client:(UIClient *)client
{
    m_client = client;
    m_viewId = [model[@"id"] integerValue];
    m_profiles = [[NSMutableDictionary alloc] init];
    m_minLength = 0;
    m_maxLength = 0;
    for (NSArray * pair in model[@"d"]) {
        NSNumber * key = pair[0];
        if (0 == m_minLength  ||  [key integerValue] < m_minLength)
            m_minLength = [key integerValue];
        if ([key integerValue] - 1 > m_maxLength)
            m_maxLength = ([key integerValue] - 1);
        Profiles * p = [[Profiles alloc] init];
        p.treated = pair[1][@"t"];
        p.untreated = pair[1][@"u"];
        [p compute];
        m_profiles[key] = p;
    }

    if (nil == m_xLabels) {
        m_xLabels = [[NSMutableArray alloc] init];
        for (NSInteger site = 0; site < m_maxLength; site += 10) {
            id<H2Label> l = [H2ViewImpl label];
            l.textColor = [H2Color blackColor];
            l.text = (0 == site ? @"S = 0" : FMT(@"%d", (int)site));
            l.alignment = NSTextAlignmentCenter;
            l.fontSize = 10;
            [self addSubview:l];
            [m_xLabels addObject:l];
        }
        m_yLabels = [[NSMutableArray alloc] init];
        for (NSInteger L = 0; L < m_maxLength; L += 10) {
            if (L < m_minLength)
                continue;
            id<H2Label> l = [H2ViewImpl label];
            l.textColor = [H2Color blackColor];
            l.text = FMT(@"%@%d", (0 == m_yLabels.count ? @"L = " : @""), (int)L);
            l.alignment = NSTextAlignmentRight;
            l.fontSize = 10;
            [self addSubview:l];
            [m_yLabels addObject:l];
        }
    }

    [self matrix_plot:model];
}

-(void)updateSiteValues:(Site *)site
{
    Profiles * p = m_profiles[@(site.L)];
    if (nil == p) {
        [site reset];
        return;
    }
    NSInteger key = site.site;
    site.treated = [p.treated[key] integerValue];
    site.untreated = [p.untreated[key] integerValue];
    site.beta = [p.beta[key] floatValue];
    site.theta = [p.theta[key] floatValue];
    site.rho = [p.rho[key] floatValue];
}

-(void)updateSite:(Site *)site withPoint:(CGPoint)point
{
    if (!CGRectContainsPoint(m_matrixFrame, point)) {
        [site reset];
        return;
    }
    CGPoint matrixLoc = CGPointMake(point.x - m_matrixFrame.origin.x, point.y - m_matrixFrame.origin.y);
    NSInteger lval = m_minLength + (matrixLoc.y / m_siteSize.height);
    NSInteger sval = matrixLoc.x / m_siteSize.width;
    if (lval < 0 || sval < 0 || sval > lval) {
        [site reset];
        return;
    }
    site.L = lval;
    site.site = sval;
    [self updateSiteValues:site];
}

-(void)layoutSubviews
{
    if (nil == m_selectedData)
        return;
    CGRect bounds = self.bounds;
    CGRect matrixFrame;
    matrixFrame.size = CGSizeMake((m_selectedData.count + m_minLength) * m_siteSize.width, m_selectedData.count * m_siteSize.height);
    matrixFrame.origin = CGPointMake((NSInteger)(0.5 * (bounds.size.width - matrixFrame.size.width)),
                                     (NSInteger)(0.5 * (bounds.size.height - matrixFrame.size.height)));
    if (!CGRectEqualToRect(matrixFrame, m_matrixFrame)) {
        m_matrixFrame = matrixFrame;
        [self startTracking:self inRect:m_matrixFrame];
    }
    H2Layout * l = [[H2Layout alloc] init];
    CGRect f = CGRectMake(0, 0, 128, 208);
    m_siteInfo.frame = [l top:f.size.height right:f.size.width of:bounds marginx:40 y:40];
    for (id<H2Label> label in m_siteLabels) {
        label.frame = [l topCentered:CGSizeMake(f.size.width - 20, 26) of:f margin:4];
        f = l.leftover;
    }
    f = m_matrixFrame;
    f.origin.y += (f.size.height + 4);
    f.size.height = 12;
    f.size.width = 26;
    f.origin.x += (0.5 * (m_siteSize.width - f.size.width));
    for (id<H2Label> label in m_xLabels) {
        label.frame = f;
        f.origin.x += 10 * m_siteSize.width;
    }
    f = m_matrixFrame;
    f.origin.x -= 32;
    f.size.width = 28;
    f.size.height = 12;
    f.origin.y += (0.5 * (m_siteSize.height - f.size.height));
    for (id<H2Label> label in m_yLabels) {
        label.frame = f;
        f.origin.y += 10 * m_siteSize.height;
    }

}

-(CGRect)frameForSite:(Site *)site
{
    CGRect siteRect = CGRectMake(m_matrixFrame.origin.x, m_matrixFrame.origin.y, m_siteSize.width, m_siteSize.height);
    siteRect.origin.y += (site.L - m_minLength) * m_siteSize.width;
    siteRect.origin.x += site.site * m_siteSize.height;
    return siteRect;
}

-(CGRect)frameForRow:(Site *)site
{
    CGRect siteRect = CGRectMake(m_matrixFrame.origin.x, m_matrixFrame.origin.y, m_siteSize.width, m_siteSize.height);
    siteRect.origin.y += (site.L - m_minLength) * m_siteSize.width;
    siteRect.size.width = (1 + site.L) * m_siteSize.width;
    return siteRect;
}

-(CGRect)frameForColumn:(Site *)site
{
    CGRect siteRect = CGRectMake(m_matrixFrame.origin.x, m_matrixFrame.origin.y, m_siteSize.width, m_siteSize.height);
    siteRect.origin.y += MAX(0, (site.site - m_minLength)) * m_siteSize.width;
    siteRect.size.height = m_matrixFrame.origin.y + (m_matrixFrame.size.height - siteRect.origin.y);
    siteRect.origin.x += site.site * m_siteSize.height;
    return siteRect;
}

-(void)handleClick:(H2Event *)event
{
    CGPoint loc = [event.context pointValue];
    [self updateSite:m_curSite withPoint:loc];
    if (!m_curSite.valid)
        return;
    H2KeyEvent * keyState = [self currentModifierKeyState];
    if (keyState.shift  ||  keyState.command) {
        [self handleSelectionClickAtSite:[m_curSite clone] withState:keyState];
        return;
    }
    if (0 < m_selection.count)
        [m_client sendMessage:@{ @"viewId" : @(m_viewId), @"event" : @"show_sel" }];
    else
        [m_client sendMessage:@{ @"viewId" : @(m_viewId), @"event" : @"site", @"L" : @(m_curSite.L), @"site" : @(m_curSite.site) }];
}

-(void)sendSelectionUpdate
{
    NSMutableArray * sites = [[NSMutableArray alloc] init];
    for (SelectionPart * part in m_selection)
        [sites addObject:@(m_plotTypeRow ? part.site.L : part.site.site)];
    [m_client sendMessage:@{ @"viewId" : @(m_viewId), @"event" : @"sel", @"data" : sites }];
}

-(void)handleSelectionClickAtSite:(Site *)site withState:(H2KeyEvent *)keyState
{
    if (nil == m_selection)
        m_selection = [[NSMutableArray alloc] init];
    SelectionPart * last = [m_selection lastObject];
    BOOL incldueIntermediates = (keyState.shift  &&  last != nil);
    NSInteger curVal = (m_plotTypeRow ? site.L : site.site);
    SelectionPart * removal = nil;
    for (SelectionPart * part in m_selection) {
        NSInteger val = (m_plotTypeRow ? part.site.L : part.site.site);
        if (val == curVal) {
            removal = part;
            break;
        }
    }

    SelectionPart * cur = removal;
    if (nil == cur)
        cur = [[SelectionPart alloc] initWithSite:site];
    NSMutableArray * partsToChange = [[NSMutableArray alloc] init];
    if (incldueIntermediates) {
        NSInteger lastVal = (m_plotTypeRow ? last.site.L : last.site.site);
        for (NSInteger s = MIN(curVal, lastVal); s <= MAX(curVal, lastVal); ++s) {
            Site * i = [site clone];
            if (m_plotTypeRow)
                i.L = s;
            else
                i.site = s;
            [partsToChange addObject:[[SelectionPart alloc] initWithSite:i]];
        }
    }
    else {
        [partsToChange addObject:cur];
    }

    if (removal) {

        for (SelectionPart * p in partsToChange) {
            NSInteger pVal = (m_plotTypeRow ? p.site.L : p.site.site);
            for (SelectionPart * c in m_selection) {
                NSInteger cVal = (m_plotTypeRow ? c.site.L : c.site.site);
                if (cVal == pVal) {
                    [c.view removeFromSuperview];
                    [m_selection removeObject:c];
                    break;
                }
            }
        }

    }
    else {

        for (SelectionPart * p in partsToChange) {
            BOOL hasAlready = NO;
            NSInteger pVal = (m_plotTypeRow ? p.site.L : p.site.site);
            for (SelectionPart * c in m_selection) {
                NSInteger cVal = (m_plotTypeRow ? c.site.L : c.site.site);
                if (cVal == pVal) {
                    hasAlready = YES;
                    break;
                }
            }
            if (!hasAlready) {
                CGRect frame;
                if (m_plotTypeRow)
                    frame = [self frameForRow:p.site];
                else
                    frame = [self frameForColumn:p.site];
                p.view.frame = frame;
                [self addSubview:p.view];
                [m_selection addObject:p];
            }
        }

    }

    [self sendSelectionUpdate];
}

-(void)resetSelection
{
    for (SelectionPart * part in m_selection)
        [part.view removeFromSuperview];
    m_selection = nil;
    [self sendSelectionUpdate];
}

typedef CGColorRef (*Color_Mapper)(CGFloat val);

Color_Mapper
get_cmap(NSString * map_name);

-(NSString *)cmapForPlotType:(NSString *)plotType
{
    static NSDictionary * s_cmapDict = nil;
    if (nil == s_cmapDict)
        s_cmapDict = @{ @"treated" : @"jet",
                        @"untreated" : @"jet",
                        @"beta" : @"viridis",
                        @"theta" : @"viridis",
                        @"rho" : @"viridis" };
    return s_cmapDict[plotType];
}

-(void)drawRect:(CGRect)dirtyRect
{
    if (0 == m_matrixFrame.origin.x)
        [self layoutSubviews];
    NSGraphicsContext * gctx = [NSGraphicsContext currentContext];
    CGContextRef ctx = (CGContextRef)[gctx graphicsPort];
    CGRect siteRect = CGRectMake(m_matrixFrame.origin.x, m_matrixFrame.origin.y, m_siteSize.width, m_siteSize.height);
    Color_Mapper cmap = get_cmap([self cmapForPlotType:m_plotType]);
    for (NSInteger L = m_minLength; ; ++L) {
        NSArray * data = m_selectedData[@(L)];
        if (nil == data)
            break;
        Profiles * prof = m_profiles[@(L)];
        for (NSInteger s = 0; s <= L; ++s) {
            CGFloat f = MIN(MAX(0.0, [data[s] floatValue] / m_max), 1.0);
            CGColorRef col = NULL;
            if (m_showNull) {
                if ([data[s] floatValue] < 0)
                    col = CGColorCreateGenericRGB(1.0, 0.0, 0.0, 1.0);
                else if (![@"treated" isEqual:m_plotType]  &&  ![@"untreated" isEqual:m_plotType]  &&
                         0 == [prof.treated[s] integerValue]  &&  0 == [prof.untreated[s] integerValue])
                    col = CGColorCreateGenericRGB(1.0, 1.0, 1.0, 1.0);
            }
            if (NULL == col)
                col = cmap(f);
            CGContextSetFillColorWithColor(ctx, col);
            CGContextFillRect(ctx, siteRect);
            CGColorRelease(col);
            siteRect.origin.x += m_siteSize.width;
        }
        siteRect.origin.x = m_matrixFrame.origin.x;
        siteRect.origin.y += m_siteSize.height;
    }
}

-(void)updateSiteInfo
{
    [m_siteLabels[0] setText:FMT(@"L: %d", (int)m_curSite.L)];
    [m_siteLabels[1] setText:FMT(@"Site: %d", (int)m_curSite.site)];
    [m_siteLabels[2] setText:FMT(@"Treated: %d", (int)m_curSite.treated)];
    [m_siteLabels[3] setText:FMT(@"Untreated: %d", (int)m_curSite.untreated)];
    [m_siteLabels[4] setText:FMT(@"beta: %0.6f", m_curSite.beta)];
    [m_siteLabels[5] setText:FMT(@"theta: %0.6f", m_curSite.theta)];
    [m_siteLabels[6] setText:FMT(@"rho: %0.6f", m_curSite.rho)];
}

-(void)mouseEnteredTrackingArea
{
}

-(void)mouseTrackingMoved:(CGPoint)point
{
    [self updateSite:m_curSite withPoint:point];
    if (!m_curSite.valid) {
        m_rowTracker.hidden = m_colTracker.hidden = YES;
        return;
    }
    m_rowTracker.frame = CGRectInset([self frameForRow:m_curSite], -1, -1);
    m_colTracker.frame = CGRectInset([self frameForColumn:m_curSite], -1, -1);
    m_rowTracker.hidden = NO;
    m_colTracker.hidden = NO;
    [self updateSiteInfo];
    m_siteInfo.hidden = NO;
}

-(void)mouseExitedTrackingArea
{
    m_siteInfo.hidden = YES;
    m_rowTracker.hidden = YES;
    m_colTracker.hidden = YES;
}

@end


@implementation Profiles

@synthesize treated = _treated;
@synthesize untreated = _untreated;
@synthesize beta = _beta;
@synthesize theta = _theta;
@synthesize rho = _rho;
@synthesize c = _c;

-(void)compute
{
    NSArray * treated = self.treated;
    NSArray * untreated = self.untreated;
    NSInteger n = treated.count - 1;
    NSMutableArray * betas = [[NSMutableArray alloc] init];
    NSMutableArray * thetas = [[NSMutableArray alloc] init];
    CGFloat treated_sum = 0.0;
    CGFloat untreated_sum = 0.0;
    CGFloat running_c_sum = 0.0;

    /* NOTE: there is an index discrepancy here between indices
       used in the code, and the indices used in the Aviran paper
       where these formulae are derived: the indices are
       reversed. so, where in the paper the formula uses
       \sum_{i=k}^{n+1}, in the code we use \sum_{i=0}^{k+1}, and
       this is intentional.

       for reference, here is the comment from the original SPATS code:
         // TargetProfile tracks an RNA of length n. arrays have n+1 entries, 
         // with index 1 corresponding to the 5'-most base, and index n 
         // corresponding to the 3'-most base.  Index 0 stores information about 
         // unmodified RNAs, where RT has fallen off the 5' end of the 
         // strand.  This convention differs from the paper, where we refer to 
         // the 5'-most base as index n, and the 3'-most base as index 1.
    */

    for (NSInteger k = 0; k < n; ++k) {
        CGFloat X_k = [treated[k] floatValue];
        CGFloat Y_k = [untreated[k] floatValue];
        treated_sum += X_k;
        untreated_sum += Y_k;
        CGFloat beta = 0.0;
        CGFloat theta = 0.0;
        if (treated_sum > 0  &&  untreated_sum > 0) {
            CGFloat Xbit = (X_k / treated_sum);
            CGFloat Ybit = (Y_k / untreated_sum);
            if (Ybit != 1.0) {
                beta = (Xbit - Ybit) / (1.0 - Ybit);
                theta = log(1.0 - Ybit) - log(1.0 - Xbit);
            }
            running_c_sum -= log(1.0 - beta);
        }
        [betas addObject:@(beta)];
        [thetas addObject:@(theta)];
    }
    [betas addObject:@(0)];
    [thetas addObject:@(0)];
    NSMutableArray * scaled_thetas = [[NSMutableArray alloc] init];
    NSMutableArray * rho = [[NSMutableArray alloc] init];
    CGFloat c = running_c_sum;
    CGFloat c_factor = (c > 0.0 ? 1.0 / c : 1.0);
    for (NSNumber * th in thetas) {
        [scaled_thetas addObject:@(c_factor * [th floatValue])];
        [rho addObject:@(((CGFloat)n) * [th floatValue])];
    }
    self.beta = betas;
    self.theta = scaled_thetas;
    self.rho = rho;
    self.c = c;
}

@end


@implementation Site
-(void)reset
{
    self.L = 0;
    self.site = -1;
    self.treated = self.untreated = 0;
    self.beta = self.theta = self.rho = 0.0;
}
-(BOOL)valid
{
    return (self.L > 0  &&  self.site >= 0);
}
-(Site *)clone
{
    Site * site = [[Site alloc] init];
    site.site = self.site;
    site.L = self.L;
    site.treated = self.treated;
    site.untreated = self.untreated;
    site.beta = self.beta;
    site.theta = self.theta;
    site.rho = self.rho;
    return site;
}
@end


@implementation SelectionPart

-(id)initWithSite:(Site *)site
{
    if ((self = [super init])) {
        _site = site;
        id<H2View> v = [H2ViewImpl container];
        v.backgroundColor = [[H2Color selectionColor] colorWithAlphaComponent:0.7];
        _view = v;
    }
    return self;
}

@end


CGFloat _viridis_data[768] = { 0.267004, 0.004874, 0.329415,
                               0.268510, 0.009605, 0.335427,
                               0.269944, 0.014625, 0.341379,
                               0.271305, 0.019942, 0.347269,
                               0.272594, 0.025563, 0.353093,
                               0.273809, 0.031497, 0.358853,
                               0.274952, 0.037752, 0.364543,
                               0.276022, 0.044167, 0.370164,
                               0.277018, 0.050344, 0.375715,
                               0.277941, 0.056324, 0.381191,
                               0.278791, 0.062145, 0.386592,
                               0.279566, 0.067836, 0.391917,
                               0.280267, 0.073417, 0.397163,
                               0.280894, 0.078907, 0.402329,
                               0.281446, 0.084320, 0.407414,
                               0.281924, 0.089666, 0.412415,
                               0.282327, 0.094955, 0.417331,
                               0.282656, 0.100196, 0.422160,
                               0.282910, 0.105393, 0.426902,
                               0.283091, 0.110553, 0.431554,
                               0.283197, 0.115680, 0.436115,
                               0.283229, 0.120777, 0.440584,
                               0.283187, 0.125848, 0.444960,
                               0.283072, 0.130895, 0.449241,
                               0.282884, 0.135920, 0.453427,
                               0.282623, 0.140926, 0.457517,
                               0.282290, 0.145912, 0.461510,
                               0.281887, 0.150881, 0.465405,
                               0.281412, 0.155834, 0.469201,
                               0.280868, 0.160771, 0.472899,
                               0.280255, 0.165693, 0.476498,
                               0.279574, 0.170599, 0.479997,
                               0.278826, 0.175490, 0.483397,
                               0.278012, 0.180367, 0.486697,
                               0.277134, 0.185228, 0.489898,
                               0.276194, 0.190074, 0.493001,
                               0.275191, 0.194905, 0.496005,
                               0.274128, 0.199721, 0.498911,
                               0.273006, 0.204520, 0.501721,
                               0.271828, 0.209303, 0.504434,
                               0.270595, 0.214069, 0.507052,
                               0.269308, 0.218818, 0.509577,
                               0.267968, 0.223549, 0.512008,
                               0.266580, 0.228262, 0.514349,
                               0.265145, 0.232956, 0.516599,
                               0.263663, 0.237631, 0.518762,
                               0.262138, 0.242286, 0.520837,
                               0.260571, 0.246922, 0.522828,
                               0.258965, 0.251537, 0.524736,
                               0.257322, 0.256130, 0.526563,
                               0.255645, 0.260703, 0.528312,
                               0.253935, 0.265254, 0.529983,
                               0.252194, 0.269783, 0.531579,
                               0.250425, 0.274290, 0.533103,
                               0.248629, 0.278775, 0.534556,
                               0.246811, 0.283237, 0.535941,
                               0.244972, 0.287675, 0.537260,
                               0.243113, 0.292092, 0.538516,
                               0.241237, 0.296485, 0.539709,
                               0.239346, 0.300855, 0.540844,
                               0.237441, 0.305202, 0.541921,
                               0.235526, 0.309527, 0.542944,
                               0.233603, 0.313828, 0.543914,
                               0.231674, 0.318106, 0.544834,
                               0.229739, 0.322361, 0.545706,
                               0.227802, 0.326594, 0.546532,
                               0.225863, 0.330805, 0.547314,
                               0.223925, 0.334994, 0.548053,
                               0.221989, 0.339161, 0.548752,
                               0.220057, 0.343307, 0.549413,
                               0.218130, 0.347432, 0.550038,
                               0.216210, 0.351535, 0.550627,
                               0.214298, 0.355619, 0.551184,
                               0.212395, 0.359683, 0.551710,
                               0.210503, 0.363727, 0.552206,
                               0.208623, 0.367752, 0.552675,
                               0.206756, 0.371758, 0.553117,
                               0.204903, 0.375746, 0.553533,
                               0.203063, 0.379716, 0.553925,
                               0.201239, 0.383670, 0.554294,
                               0.199430, 0.387607, 0.554642,
                               0.197636, 0.391528, 0.554969,
                               0.195860, 0.395433, 0.555276,
                               0.194100, 0.399323, 0.555565,
                               0.192357, 0.403199, 0.555836,
                               0.190631, 0.407061, 0.556089,
                               0.188923, 0.410910, 0.556326,
                               0.187231, 0.414746, 0.556547,
                               0.185556, 0.418570, 0.556753,
                               0.183898, 0.422383, 0.556944,
                               0.182256, 0.426184, 0.557120,
                               0.180629, 0.429975, 0.557282,
                               0.179019, 0.433756, 0.557430,
                               0.177423, 0.437527, 0.557565,
                               0.175841, 0.441290, 0.557685,
                               0.174274, 0.445044, 0.557792,
                               0.172719, 0.448791, 0.557885,
                               0.171176, 0.452530, 0.557965,
                               0.169646, 0.456262, 0.558030,
                               0.168126, 0.459988, 0.558082,
                               0.166617, 0.463708, 0.558119,
                               0.165117, 0.467423, 0.558141,
                               0.163625, 0.471133, 0.558148,
                               0.162142, 0.474838, 0.558140,
                               0.160665, 0.478540, 0.558115,
                               0.159194, 0.482237, 0.558073,
                               0.157729, 0.485932, 0.558013,
                               0.156270, 0.489624, 0.557936,
                               0.154815, 0.493313, 0.557840,
                               0.153364, 0.497000, 0.557724,
                               0.151918, 0.500685, 0.557587,
                               0.150476, 0.504369, 0.557430,
                               0.149039, 0.508051, 0.557250,
                               0.147607, 0.511733, 0.557049,
                               0.146180, 0.515413, 0.556823,
                               0.144759, 0.519093, 0.556572,
                               0.143343, 0.522773, 0.556295,
                               0.141935, 0.526453, 0.555991,
                               0.140536, 0.530132, 0.555659,
                               0.139147, 0.533812, 0.555298,
                               0.137770, 0.537492, 0.554906,
                               0.136408, 0.541173, 0.554483,
                               0.135066, 0.544853, 0.554029,
                               0.133743, 0.548535, 0.553541,
                               0.132444, 0.552216, 0.553018,
                               0.131172, 0.555899, 0.552459,
                               0.129933, 0.559582, 0.551864,
                               0.128729, 0.563265, 0.551229,
                               0.127568, 0.566949, 0.550556,
                               0.126453, 0.570633, 0.549841,
                               0.125394, 0.574318, 0.549086,
                               0.124395, 0.578002, 0.548287,
                               0.123463, 0.581687, 0.547445,
                               0.122606, 0.585371, 0.546557,
                               0.121831, 0.589055, 0.545623,
                               0.121148, 0.592739, 0.544641,
                               0.120565, 0.596422, 0.543611,
                               0.120092, 0.600104, 0.542530,
                               0.119738, 0.603785, 0.541400,
                               0.119512, 0.607464, 0.540218,
                               0.119423, 0.611141, 0.538982,
                               0.119483, 0.614817, 0.537692,
                               0.119699, 0.618490, 0.536347,
                               0.120081, 0.622161, 0.534946,
                               0.120638, 0.625828, 0.533488,
                               0.121380, 0.629492, 0.531973,
                               0.122312, 0.633153, 0.530398,
                               0.123444, 0.636809, 0.528763,
                               0.124780, 0.640461, 0.527068,
                               0.126326, 0.644107, 0.525311,
                               0.128087, 0.647749, 0.523491,
                               0.130067, 0.651384, 0.521608,
                               0.132268, 0.655014, 0.519661,
                               0.134692, 0.658636, 0.517649,
                               0.137339, 0.662252, 0.515571,
                               0.140210, 0.665859, 0.513427,
                               0.143303, 0.669459, 0.511215,
                               0.146616, 0.673050, 0.508936,
                               0.150148, 0.676631, 0.506589,
                               0.153894, 0.680203, 0.504172,
                               0.157851, 0.683765, 0.501686,
                               0.162016, 0.687316, 0.499129,
                               0.166383, 0.690856, 0.496502,
                               0.170948, 0.694384, 0.493803,
                               0.175707, 0.697900, 0.491033,
                               0.180653, 0.701402, 0.488189,
                               0.185783, 0.704891, 0.485273,
                               0.191090, 0.708366, 0.482284,
                               0.196571, 0.711827, 0.479221,
                               0.202219, 0.715272, 0.476084,
                               0.208030, 0.718701, 0.472873,
                               0.214000, 0.722114, 0.469588,
                               0.220124, 0.725509, 0.466226,
                               0.226397, 0.728888, 0.462789,
                               0.232815, 0.732247, 0.459277,
                               0.239374, 0.735588, 0.455688,
                               0.246070, 0.738910, 0.452024,
                               0.252899, 0.742211, 0.448284,
                               0.259857, 0.745492, 0.444467,
                               0.266941, 0.748751, 0.440573,
                               0.274149, 0.751988, 0.436601,
                               0.281477, 0.755203, 0.432552,
                               0.288921, 0.758394, 0.428426,
                               0.296479, 0.761561, 0.424223,
                               0.304148, 0.764704, 0.419943,
                               0.311925, 0.767822, 0.415586,
                               0.319809, 0.770914, 0.411152,
                               0.327796, 0.773980, 0.406640,
                               0.335885, 0.777018, 0.402049,
                               0.344074, 0.780029, 0.397381,
                               0.352360, 0.783011, 0.392636,
                               0.360741, 0.785964, 0.387814,
                               0.369214, 0.788888, 0.382914,
                               0.377779, 0.791781, 0.377939,
                               0.386433, 0.794644, 0.372886,
                               0.395174, 0.797475, 0.367757,
                               0.404001, 0.800275, 0.362552,
                               0.412913, 0.803041, 0.357269,
                               0.421908, 0.805774, 0.351910,
                               0.430983, 0.808473, 0.346476,
                               0.440137, 0.811138, 0.340967,
                               0.449368, 0.813768, 0.335384,
                               0.458674, 0.816363, 0.329727,
                               0.468053, 0.818921, 0.323998,
                               0.477504, 0.821444, 0.318195,
                               0.487026, 0.823929, 0.312321,
                               0.496615, 0.826376, 0.306377,
                               0.506271, 0.828786, 0.300362,
                               0.515992, 0.831158, 0.294279,
                               0.525776, 0.833491, 0.288127,
                               0.535621, 0.835785, 0.281908,
                               0.545524, 0.838039, 0.275626,
                               0.555484, 0.840254, 0.269281,
                               0.565498, 0.842430, 0.262877,
                               0.575563, 0.844566, 0.256415,
                               0.585678, 0.846661, 0.249897,
                               0.595839, 0.848717, 0.243329,
                               0.606045, 0.850733, 0.236712,
                               0.616293, 0.852709, 0.230052,
                               0.626579, 0.854645, 0.223353,
                               0.636902, 0.856542, 0.216620,
                               0.647257, 0.858400, 0.209861,
                               0.657642, 0.860219, 0.203082,
                               0.668054, 0.861999, 0.196293,
                               0.678489, 0.863742, 0.189503,
                               0.688944, 0.865448, 0.182725,
                               0.699415, 0.867117, 0.175971,
                               0.709898, 0.868751, 0.169257,
                               0.720391, 0.870350, 0.162603,
                               0.730889, 0.871916, 0.156029,
                               0.741388, 0.873449, 0.149561,
                               0.751884, 0.874951, 0.143228,
                               0.762373, 0.876424, 0.137064,
                               0.772852, 0.877868, 0.131109,
                               0.783315, 0.879285, 0.125405,
                               0.793760, 0.880678, 0.120005,
                               0.804182, 0.882046, 0.114965,
                               0.814576, 0.883393, 0.110347,
                               0.824940, 0.884720, 0.106217,
                               0.835270, 0.886029, 0.102646,
                               0.845561, 0.887322, 0.099702,
                               0.855810, 0.888601, 0.097452,
                               0.866013, 0.889868, 0.095953,
                               0.876168, 0.891125, 0.095250,
                               0.886271, 0.892374, 0.095374,
                               0.896320, 0.893616, 0.096335,
                               0.906311, 0.894855, 0.098125,
                               0.916242, 0.896091, 0.100717,
                               0.926106, 0.897330, 0.104071,
                               0.935904, 0.898570, 0.108131,
                               0.945636, 0.899815, 0.112838,
                               0.955300, 0.901065, 0.118128,
                               0.964894, 0.902323, 0.123941,
                               0.974417, 0.903590, 0.130215,
                               0.983868, 0.904867, 0.136897,
                               0.993248, 0.906157, 0.143936 };


CGColorRef
_viridis(CGFloat val)
{
    int mval = (int)((val * 255.0) + 0.5);
    return CGColorCreateGenericRGB(_viridis_data[3 * mval], _viridis_data[3 * mval + 1], _viridis_data[3 * mval + 2], 1.0);
}

CGColorRef
_grayscale(CGFloat val)
{
    return CGColorCreateGenericRGB(val, val, val, 1.0);;
}

CGFloat
_interp(CGFloat val, NSArray * d)
{
    for (int i = 0; i < d.count; ++i) {
        NSArray * e = d[i];
        if (val <= [e[0] floatValue]) {
            if (i == 0)
                return [e[1] floatValue];
            NSArray * p = d[i - 1];
            CGFloat interp = (val - [p[0] floatValue]) / ([e[0] floatValue] - [p[0] floatValue]);
            return [p[2] floatValue] + interp * ([e[1] floatValue] - [p[2] floatValue]);
        }
    }
    return [[d[d.count - 1] objectAtIndex:2] floatValue];
}

CGColorRef
_jet(CGFloat val)
{
    static NSDictionary * s_jet_data = nil;
    if (nil == s_jet_data)
        s_jet_data = @{ @"red": @[ @[ @(0.), @(0), @(0)], 
                                    @[@(0.35), @(0), @(0)],
                                    @[@(0.66), @(1), @(1)],
                                    @[@(0.89), @(1), @(1)],
                                    @[@(1), @(0.5), @(0.5)] ],
                        @"green": @[@[@(0.), @(0), @(0)],
                                   @[@(0.125), @(0), @(0)],
                                   @[@(0.375), @(1), @(1)],
                                   @[@(0.64), @(1), @(1)],
                                   @[@(0.91), @(0), @(0)],
                                   @[@(1), @(0), @(0)]],
                        @"blue":  @[@[@(0.), @(0.5), @(0.5)],
                                   @[@(0.11), @(1), @(1)],
                                   @[@(0.34), @(1), @(1)],
                                   @[@(0.65), @(0), @(0)],
                                   @[@(1), @(0), @(0)]] };
    return CGColorCreateGenericRGB(_interp(val, s_jet_data[@"red"]),
                                   _interp(val, s_jet_data[@"green"]),
                                   _interp(val, s_jet_data[@"blue"]),
                                   1.0);
}

Color_Mapper
get_cmap(NSString * map_name)
{
    if ([map_name isEqual:@"viridis"])
        return _viridis;
    else if ([map_name isEqual:@"jet"])
        return _jet;
    else
        return _grayscale;
}
