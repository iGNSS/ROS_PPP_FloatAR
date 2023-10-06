string               sat
epoch_time_msg[6]    t0          # epoch time (GPST) {eph,clk,hrclk,ura,bias,pbias}
float64[6]           udi         # SSR update interval (s) {eph,clk,hrclk,ura,bias,pbias}
uint32[6]            iod_ssr     # iod ssr {eph,clk,hrclk,ura,bias,pbias}   
uint32               iode
uint32               ura      
float64[6]           orbit_corr  # radial,along,cross,radial_dot,along_dot,cross_dot
float64[3]           clock_corr  # c0,c1,c2
float64              hrclk       # high-rate clock corection
float32[68]          code_bias
float64[68]          phase_bias

