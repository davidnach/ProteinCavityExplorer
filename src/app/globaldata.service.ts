export class GlobalData{


	pid : string;
        experimentId : string;
	
	constructor(){}


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
}







