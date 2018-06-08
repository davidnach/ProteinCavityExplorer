import { BrowserModule } from '@angular/platform-browser';
import { NgModule } from '@angular/core';
import {HttpClientModule} from '@angular/common/http';
import { FormsModule } from '@angular/forms';
import {ServerService} from './server.service';
import {UserInputService} from './userinput.service';
import {GlobalData} from './globaldata.service';
import {Routes, RouterModule} from "@angular/router";
import { ChartsModule } from 'ng2-charts';

import { Ng4LoadingSpinnerModule } from 'ng4-loading-spinner';

import { AppComponent } from './app.component';
import { UserInputComponent } from './user-input/user-input.component';
import { ExperimentResultsComponent} from './experiment-results/experiment-results.component';
import { AboutComponent } from './about/about.component';
import { HelpComponent } from './help/help.component';

import { PDBListComponent } from './PDBList/PDBList.component';

import { LoadingModule } from 'ngx-loading';
import { BarchartsComponent } from './barcharts/barcharts.component';
import {MatCheckboxModule} from '@angular/material/checkbox';
import { ChecklistModule} from 'angular-checklist';

const appRoutes: Routes = [
  {path: '',component: UserInputComponent},
  {path: 'experiment', component: ExperimentResultsComponent},
  {path: 'about', component: AboutComponent}
  {path: 'help', component: HelpComponent}
  {path: '{PDBList', component: PDBListComponent}
];

@NgModule({
  declarations: [
    AppComponent,
    UserInputComponent,
    ExperimentResultsComponent,
    AboutComponent,
HelpComponent,
    PDBListComponent,
    BarchartsComponent
  ],
  imports: [
    Ng4LoadingSpinnerModule.forRoot() ,
    LoadingModule,
    BrowserModule,
    HttpClientModule,
    FormsModule,
    RouterModule.forRoot(appRoutes),
    ChartsModule,
    MatCheckboxModule,
    ChecklistModule
  ],
  providers: [ServerService,UserInputService,GlobalData],
  bootstrap: [AppComponent]
})
export class AppModule { }
