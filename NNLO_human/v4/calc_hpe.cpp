std::array<double,25> CalcHPE_W( const Gauge& V , std::vector<int> const& idx = {} ){
    StopWatch sw , sw2;
    std::array<double,25> ans;
    std::array<double,25> time;
    SingleStaple S(V);
    mpicout << "SS: " << sw() << endl;
    sw.Reset();
    Calc_LocalElements(V,S);
    mpicout << "LP: " << sw() << endl;
    sw2.Reset();
    double sum = 0.;
    for( int i=0 ; i<25 ; i++ ){
        sw.Reset();
        const double tmp = NNLO_W( V , S , i );
	ans[i] = tmp;
	time[i] = sw();
        sum += tmp;
	if( idx.empty() )
	  mpicout << i+1 << " " << ans[i] << "  " << time[i] << " / ";
    }
    for( auto i : idx )
      mpicout << i << " " << ans[i-1] << "  " << time[i-1] << endl;
    mpicout << "total: " << sum << " " << sw2() << endl;
    sw.Reset();
    mpicout << "Total: " << NNLO_W( V , S ) << "  " << sw() << endl;
    return ans;
}

std::array<double _Complex,25> CalcHPE_P4( const Gauge& V ){
    StopWatch sw , sw2;
    std::array<double _Complex,25> ans;
    SingleStaple S(V);
    mpicout << "SS: " << sw() << endl;

    sw2.Reset();
    double _Complex sum = 0.;
    for( int i=0 ; i<25 ; i++ ){
        sw.Reset();
        const double _Complex tmp = NNLO_P4( V , S , i );
        mpicout << i+1 << " "
	    << __real(tmp) << "  " << __imag(tmp) << " "
	    << sw() << " / ";
	ans[i] = tmp;
        sum += tmp;
    }
    mpicout << endl;
    mpicout << "total: " << sum << " " << sw2() << endl;
    sw.Reset();
    mpicout << "Total: " << NNLO_P4( V , S ) << "  " << sw() << endl;
    sw.Reset();
    auto tmp = NNLO_P4( V , S , 3 )
      + NNLO_P4( V , S , 4 )
      + NNLO_P4( V , S , 5 );
    mpicout << "4+5+6: "
	    << __real(tmp) << "  " << __imag(tmp) << " "
	    << " " << sw() << endl;
    sw.Reset();
    tmp = NNLO_P4( V , S , 25 );
    mpicout << "4+5+6: "
	    << __real(tmp) << "  " << __imag(tmp) << " "
	    << " " << sw() << endl;

    sw.Reset();
    tmp = NNLO_P4( V , S , 2 ) + NNLO_P4( V , S , 19 );
    mpicout << "3+20: "
	    << __real(tmp) << "  " << __imag(tmp) << " "
	    << " " << sw() << endl;
    sw.Reset();
    tmp = NNLO_P4( V , S , 26 );
    mpicout << "3+20: "
	    << __real(tmp) << "  " << __imag(tmp) << " "
	    << " " << sw() << endl;

    sw.Reset();
    tmp = NNLO_P4( V , S , 6 ) + NNLO_P4( V , S , 23 );
    mpicout << "7+24: "
	    << __real(tmp) << "  " << __imag(tmp) << " "
	    << " " << sw() << endl;
    sw.Reset();
    tmp = NNLO_P4( V , S , 27 );
    mpicout << "7+24: "
	    << __real(tmp) << "  " << __imag(tmp) << " "
	    << " " << sw() << endl;

    tmp = NNLO_P4( V , S , 20 )
      + NNLO_P4( V , S , 21 )
      + NNLO_P4( V , S , 22 );
    mpicout << "21+22+23: "
	    << __real(tmp) << "  " << __imag(tmp) << " "
	    << " " << sw() << endl;
    sw.Reset();
    tmp = NNLO_P4( V , S , 3 )
      + NNLO_P4( V , S , 4 )
      + NNLO_P4( V , S , 5 )
      + NNLO_P4( V , S , 20 )
      + NNLO_P4( V , S , 21 )
      + NNLO_P4( V , S , 22 );
    mpicout << "4+5+6+21+22+23: "
	    << __real(tmp) << "  " << __imag(tmp) << " "
	    << " " << sw() << endl;
    sw.Reset();
    tmp = NNLO_P4( V , S , 28 );
    mpicout << "4+5+6+21+22+23: "
	    << __real(tmp) << "  " << __imag(tmp) << " "
	    << " " << sw() << endl;

    return ans;
}
