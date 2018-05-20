export class GlobalData{


	pid : string;
  experimentId : string;
	pocketToResidues;

	constructor(){
		this.pocketToResidues = new Map<number,number[]>();
	}


	setPid(pid : string){
		this.pid = pid;
	}

	getPid(){
		return this.pid;
	}

	setExperimentId(experimentId : string){
		this.experimentId = experimentId;
	}

	getExperimentId(){
		return this.experimentId;

	}
	
	getPocketNums(){
	       return this.pocketToResidues.keys();
	}

	getMap(){
		return this.pocketToResidues;
	}

	getPocketResidues(pocketNum : number){
		return this.pocketToResidues.get(pocketNum);
	}

	setPockets(pockets){
		var residuesString;
		var residues = [];
		var i;
		var j;
		for(i = 0; i < pockets.length; i++){
			residuesString = pockets[i].residues.split(',');
			for (j = 0; j < residuesString.length; j++){
				residues.push(Number(residuesString[j]));
			}
			this.pocketToResidues.set(i,residues);
			residues = [];
		}
		console.log(this.pocketToResidues);
	}
}
