
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
@end


@interface SpatsMatrix : H2ViewImpl_View < H2ViewDrawer >
{
    UIClient * m_client;
    NSInteger m_viewId;
    NSMutableDictionary * m_profiles;
    NSInteger m_minLength;
    CGSize m_siteSize;
    CGRect m_matrixFrame;
    Site * m_curSite;

    NSMutableDictionary * m_selectedData;
    CGFloat m_max;
}

@end


@implementation SpatsMatrix

-(id)init
{
    if ((self = [super init])) {
        self.drawer = self;
        m_curSite = [[Site alloc] init];
        m_siteSize = CGSizeMake(5, 5);
        m_matrixFrame = CGRectMake(0, 0, 800, 800);
        [self addClickHandler:[H2EventHandler handlerWithTarget:self selector:@selector(handleClick:)]];
    }
    return self;
}

-(void)matrix_plot:(NSDictionary *)params
{
    m_max = [params[@"max"] floatValue];
    m_selectedData = [[NSMutableDictionary alloc] init];
    SEL prof_sel = NSSelectorFromString(params[@"plot"]);
    for (id key in m_profiles.allKeys) {
        Profiles * p = m_profiles[key];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Warc-performSelector-leaks"
        m_selectedData[key] = [p performSelector:prof_sel];
#pragma clang diagnostic pop
    }
    [self setNeedsDisplay];
}

-(void)updateWithModel:(NSDictionary *)model client:(UIClient *)client
{
    m_client = client;
    m_viewId = [model[@"id"] integerValue];
    m_profiles = [[NSMutableDictionary alloc] init];
    m_minLength = 0;
    for (NSArray * pair in model[@"d"]) {
        NSNumber * key = pair[0];
        if (0 == m_minLength  ||  [key integerValue] < m_minLength)
            m_minLength = [key integerValue];
        Profiles * p = [[Profiles alloc] init];
        p.treated = pair[1][@"t"];
        p.untreated = pair[1][@"u"];
        [p compute];
        m_profiles[key] = p;
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

-(void)handleClick:(H2Event *)event
{
    CGPoint loc = [event.context pointValue];
    [self updateSite:m_curSite withPoint:loc];
    if (!m_curSite.valid)
        return;
    [m_client sendMessage:@{ @"viewId" : @(m_viewId), @"event" : @"site", @"L" : @(m_curSite.L), @"site" : @(m_curSite.site) }];
}

-(void)drawRect:(CGRect)dirtyRect
{
    CGRect bounds = self.bounds;
    m_matrixFrame.size = CGSizeMake((m_selectedData.count + m_minLength) * m_siteSize.width, m_selectedData.count * m_siteSize.height);
    m_matrixFrame.origin = CGPointMake((NSInteger)(0.5 * (bounds.size.width - m_matrixFrame.size.width)),
                                       (NSInteger)(0.5 * (bounds.size.height - m_matrixFrame.size.height)));
    NSGraphicsContext * gctx = [NSGraphicsContext currentContext];
    CGContextRef ctx = (CGContextRef)[gctx graphicsPort];
    CGRect siteRect = CGRectMake(m_matrixFrame.origin.x, m_matrixFrame.origin.y, m_siteSize.width, m_siteSize.height);
    for (NSInteger L = m_minLength; ; ++L) {
        NSArray * data = m_selectedData[@(L)];
        if (nil == data)
            break;
        for (NSInteger s = 0; s <= L; ++s) {
            CGFloat f = MIN(MAX(0.0, [data[s] floatValue] / m_max), 1.0);
            CGColorRef col = CGColorCreateGenericRGB(f, f, f, 1.0);
            CGContextSetFillColorWithColor(ctx, col);
            CGContextFillRect(ctx, siteRect);
            CGColorRelease(col);
            siteRect.origin.x += m_siteSize.width;
        }
        siteRect.origin.x = m_matrixFrame.origin.x;
        siteRect.origin.y += m_siteSize.height;
    }
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
                beta = MAX(0.0, (Xbit - Ybit) / (1.0 - Ybit));
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
    self.beta = self.theta = self.rho = 0.0;
}
-(BOOL)valid
{
    return (self.L > 0  &&  self.site >= 0);
}
@end
