class VertexPair {
    constructor(a,b){
        this.v1 = a;
        this.v2 = b;
    }
    static equals(e1,e2){
        if(e1.v1==e2.v1 && e1.v2==e2.v2){
            return true;
        }
        if(e1.v1==e2.v2 && e1.v2==e2.v1){
            return true;
        }
        return false;
    }
}

s = [];

vp1 = new VertexPair(1,2);
vp2 = new VertexPair(1,2);
vp3 = new VertexPair(3,2);

temp = [vp1, vp2, vp3];

for(var i=0;i<temp.length;i++){
    console.log(temp[i]);
    if(s.length==0){
        s.push(temp[i]);
    }
    else{
        for(var j=0;j<s.length;j++){
            console.log(s[j]);
            if(VertexPair.equals(s[j], temp[i])){
                break;
            }
            s.push(temp[i]);
        }
    }
    
}

console.log(s);