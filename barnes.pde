/* #defines */
int PARTICLES       = 1500; // how many stars
int PARTICLE_RADIUS = 2;    // how big are they
int GRID_SIZE       = 800;
float G             = 6.67e-10;  // gravitational constant (adjusted)

float MASS_UPPER_LIMIT  = 2.5; // what's the largest mass we'll allow
float MASS_SCALE_FACTOR = 1e4;  // 1 solar mass
float SEP_SCALE_FACTOR  = 1e2; // 1 AU
float TIME_SCALE_FACTOR = 5e3;  // how long is one timestep?
float SOFT_COEFFICIENT  = 25*SEP_SCALE_FACTOR; 
float BLACK_HOLE_MASS   = 4.31e6*MASS_SCALE_FACTOR; // from wikipedia
float INIT_V= 100.0;      // initial velocity scaling factor for stars

float YSKEW = 0.15;      // Smaller skew will compress the galaxy vertically
int   USE_EXP = 1;       // use an exponentially distributed galaxy
float LAMBDA_INIT = 3.5; // effects lambda of the exponential distribution

ArrayList star_colors;
                 
class Rgb {
  int r;
  int g;
  int b;
  char sp5;
  
  Rgb (int _r, int _g, int _b, char _s) {
    r = _r;
    b = _b;
    g = _g;
    sp5 = _s;
  }
}
  

float sample_exp (float a, float lambda) {
  float u = random(1.0);
  return a*(log(u)/-lambda);
}


PBody[] qs = new PBody[PARTICLES+1];

class PVector {
    float x;
    float y;
    
    PVector (float _x, float _y) {
        x = _x;
        y = _y;
    }
    
    void sub (PVector v) {
        x -= v.x;
        y -= v.y;
    }
    
    void add (PVector v) {
        x += v.x;
        y += v.y;
    }

    void scalar_mult (float _a) {
        x *= _a;
        y *= _a;
    }
    
    float dot (PVector v) {
        return v.y*y + v.x*x;
    }
    
    float mag() {
        return sqrt((float)(x*x) + (float)(y*y));
    }
    
}

class PBody {

    PVector v;
    PVector pos;
    float m;
    float rad;
    int idx;
    Rgb star_type;
    
    PBody(float vx, float vy, float posx, float posy, float _m, float _r, int _i, Rgb _s) {
        v = new PVector(vx, vy);
        pos = new PVector(posx, posy);
        m = _m;
        rad = _r;
        idx = _i;
        star_type = _s;
    }

    PVector get_pos() {
        return pos;
    }

    void set_pos(PVector _p) {
        pos = _p;
    }

    PVector get_v() {
        return v;
    }

    void set_v(PVector _v) {
        v = _v;
    }

    float get_mass() {
        return m;
    }
    
    float get_rad() {
        return rad;
    }
    
    void set_rad(float _r) {
        rad = _r;
    }

    void set_mass(float _m) {
        m = _m;
    }

    void update() {
        PVector f = new PVector(0, 0);
      
        for (int i = 0; i < PARTICLES+1; i++) {
          
            if (idx == i) {
              continue;
            }

            PVector r = new PVector((qs[i].pos.x-pos.x)*SEP_SCALE_FACTOR,
                                    (qs[i].pos.y-pos.y)*SEP_SCALE_FACTOR);
            

          
            float fx;
            float fy;
            float rmag = r.mag();
            float mass_o = qs[i].get_mass();
            float rmag2 = rmag * rmag;
            float sum = rmag2 + SOFT_COEFFICIENT * SOFT_COEFFICIENT;
            float denom = sum * sqrt(sum);
            
              //float rmag3 = rmag * rmag * rmag;
             fx = (G * m * mass_o * r.x) / denom;
             fy = (G * m * mass_o * r.y) / denom;
              
 
            f.add(new PVector(fx, fy));
        }

        // F = ma => a = F/m
        // v = at = F/m * t 
        // t == 1 ==> v_n = v_i + (F/m)*t

        v.add(new PVector(TIME_SCALE_FACTOR*(f.x/m), TIME_SCALE_FACTOR*(f.y/m)));
      
        PVector nv = new PVector(TIME_SCALE_FACTOR/SEP_SCALE_FACTOR*v.x, TIME_SCALE_FACTOR/SEP_SCALE_FACTOR*v.y);
        pos.add(nv);

    }
    
    void draw() {
        if (idx != PARTICLES) {
          noStroke();
          float _x = pos.x - width/2;
          float _y = pos.y - height/2;
          float x = sqrt(width/2*width/2 + height/2*height/2)/sqrt(_x*_x + _y*_y);
          fill(star_type.r,star_type.g,star_type.g, 255*x);
          ellipse((float)pos.x, (float)pos.y, (float)rad, (float)rad);
        } else {
          noStroke();
          fill(160, 43, 219, 140);
          ellipse((float)pos.x, (float)pos.y, 4, 4);
        }
    }
}


void setup() {

    size(800, 800);
    background(0);
    float a = 0.0;
    float inc = TWO_PI/80.0;
    int i;
    
    star_colors = new ArrayList();
star_colors.add(new Rgb(155,176,255,'O'));
star_colors.add(new Rgb(170,191,255,'B'));
star_colors.add(new Rgb(202,215,255,'A'));
star_colors.add(new Rgb(248,247,255,'F'));
star_colors.add(new Rgb(255,244,234,'G'));
star_colors.add(new Rgb(255,210,151,'K'));
star_colors.add(new Rgb(255,204,111,'M'));

    for (i = 0; i < PARTICLES; i++) {
        float mass = random(1.0, MASS_UPPER_LIMIT) * MASS_SCALE_FACTOR;
        float radius = PARTICLE_RADIUS * (mass/MASS_SCALE_FACTOR);
        float theta = random(TWO_PI);
        float placement;
        if (USE_EXP == 1) {
          placement = sample_exp(1.0, LAMBDA_INIT);
        } else {
          placement = random(0.25);
        }
        float galposx = (width/2) + (width/2)*placement*cos(theta);
        float galposy = (height/2) + (height/2)*YSKEW*placement*sin(theta);
        float difx = (galposx - width/2)*SEP_SCALE_FACTOR;
        float dify = (galposy - height/2)*SEP_SCALE_FACTOR;
        float galvx = -INIT_V*(1/sqrt(difx*difx + dify*dify))*sin(atan2(dify, difx));
        float galvy = INIT_V*(1/sqrt(difx*difx + dify*dify))*cos(atan2(dify,difx));
       
        qs[i] = new PBody(galvx,  // vx
                          galvy,  // vy
                          galposx,
                          galposy, 
                          mass,
                          radius,
                          i,
                          (Rgb)star_colors.get(6 - (int)sample_exp(6.0, 1.0)%7));
                  /* 
         qs[i] = new PBody(0,  // vx
                          0,  // vy
                          random(width/4, 3*width/4), 
                          random(height/4, 3*height/4),
                          mass,
                          radius,
                          i);
                          */
        a += inc;
    }

    /* add a black hole for the mouse */
    qs[i] = new PBody(0, 0, 
              width/2, height/2,
              BLACK_HOLE_MASS, 
              29511.0,
              i,
              new Rgb(255,255,255,'z'));

}

void mouseMoved() {
    PBody bh = qs[PARTICLES];
    bh.pos.x = mouseX;
    bh.pos.y = mouseY;
}

void draw() {
  
    background(0);
    
    for (int i = 0; i < PARTICLES+1; i++) {
        qs[i].update();
        qs[i].draw();
    }
}