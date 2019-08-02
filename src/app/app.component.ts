import { Component } from '@angular/core';
import {Router} from '@angular/router';
import {ServerService} from './server.service';
import {UserInputService} from './userinput.service';
import {GlobalData} from './globaldata.service';

@Component({
  selector: 'app-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.css'],
  providers :[GlobalData]
})

export class AppComponent {
     public loading = false;  
     constructor(private serverService: ServerService,
     private router: Router,private userInputService: UserInputService,
     private globaldata : GlobalData ){
     
     this.userInputService.pidEntered.subscribe(
        (pid : string) => {
            this.submitPidToServer(pid);
     });
     this.userInputService.experimentIdEntered.subscribe(
	      (expId : string) => {
	          this.submitExperimentId(expId);
	   });
 }

  submitPidToServer(pid : string){
      this.loading = true;
      this.serverService.createExperiment({pid: pid})
          .subscribe((response) => {
              console.log(response['expId']);
              this.globaldata.setExperimentId(response['expId']);
              this.globaldata.setPockets(response['pockets']);
              this.loading = false;
              console.log(this.loading);	
              this.router.navigate(['experiment']);},
              (error) => {
                  this.loading = false;
                  console.log(error);
               }
        );
  }
  
  submitExperimentId(experimentId : string){
      this.loading = true;
      this.serverService.retrieveExperiment({expId : experimentId})
      .subscribe((response) => {
          console.log(response);
          if(response === 'invalid') {
              this.globaldata.setExperimentId('invalid');
          } else {
              this.globaldata.setPid(response['pid']);
              this.globaldata.setExperimentId(response['expId']);
              this.globaldata.setPockets(response['pockets']);
              this.loading = false;
              this.router.navigate(['experiment']);
          }	
      },	
      (error) => {
          console.log(error)
          this.loading = false;
      }
   );
 }

}
