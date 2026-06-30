std::array<double,25> CalcHPE_W( const Gauge& V ){
    StopWatch sw , sw2;
    std::array<double,25> ans;
    SingleStaple S(V);
    mpicout << "SS: " << sw() << endl;
    sw.Reset();
    Calc_LocalPlaquette(V,S);
    mpicout << "LP: " << sw() << endl;
    sw2.Reset();
    double sum = 0.;
    for( int i=0 ; i<25 ; i++ ){
        sw.Reset();
        const double tmp = NNLO_W( V , S , i );
        mpicout << i+1 << " " << tmp << "  " << sw() << " / ";
	ans[i] = tmp;
        sum += tmp;
    }
    mpicout << "total: " << sum << " " << sw2() << endl;
    sw.Reset();
    mpicout << "Total: " << NNLO_Wall( V , S ) << "  " << sw() << endl;
    mpicout << "end" << endl;
    return ans;
}

std::array<double,25> CalcHPE_P4( const Gauge& V ){
    StopWatch sw , sw2;
    std::array<double,25> ans;
    SingleStaple S(V);
    mpicout << "SS: " << sw() << endl;

    sw2.Reset();
    double sum = 0.;
    for( int i=0 ; i<25 ; i++ ){
        sw.Reset();
        const double tmp = NNLO_P4( V , S , i );
        mpicout << i+1 << " " << tmp << "  " << sw() << " / ";
	ans[i] = tmp;
        sum += tmp;
    }
    mpicout << "total: " << sum << " " << sw2() << endl;
    sw.Reset();
    mpicout << "Total: " << NNLO_P4all( V , S ) << "  " << sw() << endl;
    mpicout << "end" << endl;
    return ans;
}
