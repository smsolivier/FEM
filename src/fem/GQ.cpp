#include "GQ.H"
#include <cassert>
#include <iostream>

GQ::GQ(int a_p) {

	m_p = a_p; 

	if (a_p == 1) {
		m_w = {
			1.0
		}; 

		m_xy = {
			{0.33333333333333333333,  0.33333333333333333333}
		}; 
	}

	else if (a_p == 2) {
		m_w = {
			0.33333333333333333333, 
			0.33333333333333333333, 
			0.33333333333333333333 
		}; 

		m_xy = {
			{0.0,  0.5}, 
			{0.5,  0.0}, 
			{0.5,  0.5}
		}; 
	}

	else if (a_p == 3) {
		m_w = {
			0.30000000000000000000, 
			0.30000000000000000000, 
			0.30000000000000000000, 
			0.033333333333333333333, 
			0.033333333333333333333, 
			0.033333333333333333333
		}; 

		m_xy = {
			{0.66666666666666666667,  0.16666666666666666667}, 
			{0.16666666666666666667,  0.66666666666666666667}, 
			{0.16666666666666666667,  0.16666666666666666667}, 
			{0.0,  0.5}, 
			{0.5,  0.0}, 
			{0.5,  0.5} 
		}; 
	}

	else if (a_p == 4) {
		m_w = {
			0.22338158967801146570, 
			0.22338158967801146570, 
			0.22338158967801146570, 
			0.10995174365532186764, 
			0.10995174365532186764, 
			0.10995174365532186764
		}; 

		m_xy = {
			{0.10810301816807022736,  0.44594849091596488632}, 
			{0.44594849091596488632,  0.10810301816807022736}, 
			{0.44594849091596488632,  0.44594849091596488632}, 
			{0.81684757298045851308,  0.091576213509770743460}, 
			{0.091576213509770743460,  0.81684757298045851308}, 
			{0.091576213509770743460,  0.091576213509770743460}
		}; 
	}

	else if (a_p == 5) {
		m_w = {
			0.12593918054482715260, 
			0.12593918054482715260, 
			0.12593918054482715260, 
			0.13239415278850618074, 
			0.13239415278850618074, 
			0.13239415278850618074, 
			0.22500000000000000000
		}; 

		m_xy = {
			{0.79742698535308732240,  0.10128650732345633880}, 
			{0.10128650732345633880,  0.79742698535308732240}, 
			{0.10128650732345633880,  0.10128650732345633880}, 
			{0.059715871789769820459,  0.47014206410511508977}, 
			{0.47014206410511508977,  0.059715871789769820459}, 
			{0.47014206410511508977,  0.47014206410511508977}, 
			{0.33333333333333333333,  0.33333333333333333333}
		}; 
	}

	else if (a_p == 6) {
		m_w = {
			0.050844906370206816921,
			0.050844906370206816921,
			0.050844906370206816921,
			0.11678627572637936603,
			0.11678627572637936603,
			0.11678627572637936603,
			0.082851075618373575194,
			0.082851075618373575194,
			0.082851075618373575194,
			0.082851075618373575194,
			0.082851075618373575194,
			0.082851075618373575194 
		};

		m_xy = {
			{0.87382197101699554332,   0.063089014491502228340}, 
			{0.063089014491502228340,  0.873821971016995543320}, 
			{0.063089014491502228340,  0.063089014491502228340}, 
			{0.50142650965817915742,   0.249286745170910421290}, 
			{0.24928674517091042129,   0.501426509658179157420}, 
			{0.24928674517091042129,   0.249286745170910421290}, 
			{0.053145049844816947353,  0.310352451033784405420}, 
			{0.31035245103378440542,   0.053145049844816947353}, 
			{0.053145049844816947353,  0.636502499121398647230}, 
			{0.31035245103378440542,   0.636502499121398647230}, 
			{0.63650249912139864723,   0.053145049844816947353}, 
			{0.63650249912139864723,   0.310352451033784405420}
		}; 
	}

	else if (a_p == 8) {
		m_w = {
			0.14431560767779,
			0.09509163426728,
			0.09509163426728,
			0.09509163426728,
			0.10321737053472,
			0.10321737053472,
			0.10321737053472,
			0.03245849762320,
			0.03245849762320,
			0.03245849762320,
			0.02723031417443,
			0.02723031417443,
			0.02723031417443,
			0.02723031417443,
			0.02723031417443,
			0.02723031417443
		}; 

		m_xy = {
			{0.33333333333333, 0.33333333333333}, 
			{0.45929258829272, 0.45929258829272}, 
			{0.45929258829272, 0.08141482341455}, 
			{0.08141482341455, 0.45929258829272}, 
			{0.17056930775176, 0.17056930775176}, 
			{0.17056930775176, 0.65886138449648}, 
			{0.65886138449648, 0.17056930775176}, 
			{0.05054722831703, 0.05054722831703}, 
			{0.05054722831703, 0.89890554336594}, 
			{0.89890554336594, 0.05054722831703}, 
			{0.26311282963464, 0.72849239295540}, 
			{0.72849239295540, 0.00839477740996}, 
			{0.00839477740996, 0.26311282963464}, 
			{0.72849239295540, 0.26311282963464}, 
			{0.26311282963464, 0.00839477740996}, 
			{0.00839477740996, 0.72849239295540}
		}; 
	}

	else {
		cout << "ERROR: GQ order not defined" << endl; 
		exit(0); 
	}

	assert(m_w.size() == m_xy.size()); 
	m_Np = m_w.size(); 
}